import sys
import os

try:
    sys.argv[1]
except IndexError:
    datafile = 'data.lammps'
else:
    datafile = sys.argv[1]

try:
    sys.argv[2]
except IndexError:
    linkfile = 'polyout'
else:
    linkfile = sys.argv[2]

#print "'python cross_link_density.py ",datafile
#os.system('python cross_link_density.py '+datafile) 

n_tetra = 3000
n_amine = n_tetra

a1 = 0
a2 = 0

with open(linkfile,'r') as f:
    for line in f:
        line = line.split()
        if line[0] in ['7','14'] and line[1] in ['7','14']:
            a1 += int(line[4])
        if line[0] in ['7','22'] and line[1] in ['7','22']:
            a2 += int(line[4])

print 'epoxy amine1:',a1,' / ',n_amine*2 ,'\t',100* float(a1)/(n_amine*2),'%'
print 'epoxy amine2:',a2,' / ',n_amine*2,'\t',100* float(a2)/(n_amine*2),'%'
print 100 * (a1+a2) / (n_amine*4), '% cross linked'
