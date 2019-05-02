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
vdw_defs = str(sys.argv[3])

lammps_args = ['aprun', '-n', str(n_cores), 'lmp_xc30', '-in'] # default archer mpi archer args

lammps_dir = 'lammps_inputs/'

primary_amine      = json.load(open('configs/epoxy_amine1.json'))
secondary_amine    = json.load(open('configs/epoxy_amine2.json'))
close_unused_epoxy = json.load(open('configs/epoxy_close.json'))

compress = True # do you want to run a lammps compression before crosslinking

n_possible_bonds = 400 # if 100% crosslinking was achieved
crosslink_density = 0.80 # desired crosslink density
n_desired_bonds  = int(n_possible_bonds * crosslink_density) # desired number of crosslink bonds

bonds_per_loop = 20
search_radius = 6
mask_radius = 20

outputfile = 'polyout'
networkfile = 'nx.dat'

def react(sim, reaction_config, n_attempt_bonds, new_connections):
  makeBond = MakeBond(sim, reaction_config, 
                Nbonds = n_attempt_bonds,
                new_connections=new_connections,
                search_radius = search_radius,
                outputfile = outputfile,
                networkfile = networkfile,
                mask_radius = mask_radius)
  return makeBond.bonds_made, makeBond.new_connections
  
def timer(func):
    start = time.time()
    def wrapper(*args,**kwargs):
        return func(*args,**kwargs)
    end = time.time()
    print 'Total wall time: ',end-start,'seconds'
    return wrapper

@timer
def polymerise(reaction_config, n_attempt_bonds, filename):
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
  return n_new_bonds

def loop_polymerise(reaction, bonds_made):
  # crosslinks a file called 'to_be_crosslinked.data' as much as possible according to 
  # specifiaction; outputs 'reacted.data'
  assert(os.path.isfile('to_be_crosslinked.data'))

  while bonds_made < n_desired_bonds:
    diff = n_desired_bonds - bonds_made
    if diff < bonds_per_loop:
      n_attempt_bonds = diff
    else:
      n_attempt_bonds = bonds_per_loop

    # Makes some crosslins on the to_be_crosslinked.data simulation, outputs crosslinked.data
    n_new_bonds = polymerise(reaction, n_attempt_bonds, 'to_be_crosslinked.data')
    bonds_made += n_new_bonds
    
    if n_new_bonds:
      # Runs a lammps minimization on crosslinked.data, outputs minimized.data
      return_code = subprocess.call(lammps_args + [lammps_dir+'in.minimize'])
      if return_code:
        print 'Compress failed';  sys.exit()
    else:
      # structure unchanged 
      subprocess.call(['mv','crosslinked.data','minimized.data'])

    if (   (bonds_made == n_desired_bonds)
        or (n_new_bonds == 0)):
      subprocess.call(['mv','minimized.data','reacted.data'])      
      break
    else:
      # more polymerisations to be done in this loop
      subprocess.call(['mv','minimized.data','to_be_crosslinked.data'])      
      
  return bonds_made

# Main control sequence

if compress:
  # compress monomers to atmospheric pressure
  subprocess.call(['cp',inputfile,'box.data']) # so in.compress reads the correct file
  return_code = subprocess.call(lammps_args + [lammps_dir+'in.compress'])
  if return_code:
    print 'Compress failed';  sys.exit()
  subprocess.call(['cp','compressed.data','equilibrated.data']) # prepare lammps file for loop
else:
  subprocess.call(['cp',inputfile,'equilibrated.data'])

bonds_made = 0
loop = 0
while bonds_made < n_desired_bonds:
  loop += 1
  
  subprocess.call(['mv','equilibrated.data','to_be_crosslinked.data'])
  bonds_made = loop_polymerise(primary_amine, bonds_made)

  subprocess.call(['mv','reacted.data','to_be_crosslinked.data'])
  bonds_made = loop_polymerise(secondary_amine, bonds_made)

  # Runs a lammps npt on reacted.data, outputs new_box.data
  return_code = subprocess.call(lammps_args + [lammps_dir+'in.equil'])
  if return_code:
    print 'Equilibration failed',loop;  sys.exit()
  
mask_radius = 0 # all unused epoxys need to closed at once
polymerise(close_unused_epoxy, n_possible_bonds, 'equilibrated.data')

return_code = subprocess.call(lammps_args + [lammps_dir+'in.minimize'])
if return_code: print 'Equilibration failed',loop;  sys.exit()

subprocess.call(lammps_args + [lammps_dir+'in.long_equil'])
if return_code: print 'Equilibration failed',loop;  sys.exit()
    
subprocess.call(['mv','equilibrated.data','final.data'])
