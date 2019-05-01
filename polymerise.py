import subprocess
import time
import sys, os
import json
from scripts import ReadLammpsData
from scripts import MakeBond
from scripts import Writer
from forcefield import Parameterise

""" 
Script to manage an itterative python crosslinking algorithm and lammps simulation
Input file is compressed to atmospheric conditions
Primary and Secondary amines are reacted in that order untill desired bonds made

python polymerisation.py n_cores inputfile
"""

n_cores = int(sys.argv[1])
inputfile = str(sys.argv[2])
#n_cores = '192'
lammps_args = ['aprun', '-n', str(n_cores), 'lmp_xc30', '-in']
lammps_dir = 'lammps_inputs/'
primary_amine      = json.load(open('configs/epoxy_amine1.json'))
secondary_amine    = json.load(open('configs/epoxy_amine2.json'))
close_unused_epoxy = json.load(open('configs/epoxy_close.json'))

n_possible_bonds = 400 # if 100% crosslinking was achieved
n_desired_bonds  = 320  # desired crosslink bonds before end

bonds_per_loop = 20
search_radius = 6

outputfile = 'polyout'
networkfile = 'nx.dat'
vdw_defs = 'vdw.json'

def react(sim, reaction_config, n_attempt_bonds, new_connections):
  makeBond = MakeBond(sim, reaction_config, 
                Nbonds = n_attempt_bonds,
                new_connections=new_connections,
                search_radius = search_radius,
                outputfile = outputfile,
                networkfile = networkfile)
  return makeBond.bonds_made, makeBond.new_connections

def polymerise(reaction_config, n_attempt_bonds, filename):
  start = time.time() # only needed for benchmarking

  sim = ReadLammpsData(filename)

  # If this is part of an itterative polymerisation, check if new connecions
  # have been made so they don't need to be parameterised every time
  if os.path.isfile('new_connections.json'):
    new_connections = 'new_connections.json'
  else:
    new_connections = {'angles':{},'dihedrals':{}}

  sim.vdw_defs = json.load(open( vdw_defs ))
  sim.vdw_defs = {int(k):v for k,v in sim.vdw_defs.items()}

  n_new_bonds, new_connections = react(sim, reaction_config, n_attempt_bonds, new_connections)
  
  Parameterise(sim, sim.vdw_defs, new_connections)

  json.dump(sim.vdw_defs,open('vdw.json','w'))
  output = Writer(sim)
  output.write_lammps('crosslinked.data')
  end = time.time()
  print 'Total wall Time for crosslinking',end - start
  return n_new_bonds

def loop_polymerise(reaction, bonds_made):
  while bonds_made < n_desired_bonds:
    diff = n_desired_bonds - bonds_made
    if diff < bonds_per_loop:
      n_attempt_bonds = diff
    else:
      n_attempt_bonds = bonds_per_loop

    # Makes some crosslins on the minimized.data simulation, outputs crosslinked.data
    n_new_bonds = polymerise(reaction, n_attempt_bonds, 'minimized.data')
    bonds_made += n_new_bonds
    subprocess.call(['rm','minimized.data']) # just to be safe

    # Runs a lammps minimization on crosslinked.data, outputs minimized.data
    return_code = subprocess.call(lammps_args + [lammps_dir+'in.minimize'])
    if return_code:
      print 'Compress failed';  sys.exit()
    subprocess.call(['rm','crosslinked.data']) # just to be safe

    if (   (bonds_made == n_desired_bonds)
        or (n_new_bonds == 0)):
      break

  return bonds_made



# compress random arrangement of monomers to atmospheric pressure
subprocess.call(['cp',inputfile,'box.data']) # so in.compress reads the correct file
return_code = subprocess.call(lammps_args + [lammps_dir+'in.compress'])
if return_code:
  print 'Compress failed';  sys.exit()
subprocess.call(['cp','compressed.data','minimized.data']) # prepare lammps file for loop

#subprocess.call(['cp',inputfile,'minimized.data']) # for testing

bonds_made = 0
loop = 0
while bonds_made < n_desired_bonds:
  loop += 1
  bonds_made = loop_polymerise(primary_amine, bonds_made)
  bonds_made = loop_polymerise(secondary_amine, bonds_made)

  # Runs a lammps npt on minimized.data, outputs new_box.data
  if loop < 5:
    return_code = subprocess.call(lammps_args + [lammps_dir+'in.equil'])
  else:
    return_code = subprocess.call(lammps_args + [lammps_dir+'in.hot_equil'])
  if return_code:
    print 'Compress failed';  sys.exit()
  
  subprocess.call(['cp','new_box.data',str(bonds_made)+'.data']) 
  
polymerise(close_unused_epoxy, n_possible_bonds, 'new_box.data')
subprocess.call(lammps_args + [lammps_dir+'in.long_equil'])
    

