import subprocess
import time
import sys, os
import json
from scripts import ReadLammpsData
from scripts import MakeBond
from scripts import Writer
from forcefield import Parameterise

from Cython.Build import cythonize


""" 
Script to manage an itterative python crosslinking algorithm and lammps simulation
Input file is compressed to atmospheric conditions
Primary and Secondary amines are reacted in that order untill desired bonds made

python polymerisation.py n_cores inputfile
"""

vdw_defs = 'vdw.json'


lammps_dir = 'lammps_inputs/'

primary_amine      = json.load(open('configs/epoxy_amine1.json'))
secondary_amine    = json.load(open('configs/epoxy_amine2.json'))
close_unused_epoxy = json.load(open('configs/epoxy_close.json'))

compress = True # do you want to run a lammps compression before crosslinking

crosslink_density = 0.80 # desired crosslink density

bonds_per_loop = 160
search_radius = 6
mask_radius = 0

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
 

polymerise(primary_amine , 200, sys.argv[1])

