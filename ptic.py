#!/usr/bin/env python3
# ptic.py
# A script to compute the p-tic order parameter based on a given neighbour
# list. Here, the local (i.e., per point) p-tic order parameter is defined as
# follows:
#
#   psi_{p,j} = \frac{1}{N_{j,nn}}\sum_{k\in nn}\exp(ip\theta_{jk})
#
# where N_{j,nn} is the number of nearest neighbours for point j and
# \theta_{jk} is the angle that the bond vector connecting points j and k
# makes with the x-axis.

import sys
import numpy as np

args = sys.argv
if (len(args) != 4):
    print("Usage: ptic.py p neigh_file out_file")
    sys.exit(0)

p = int(args.pop(1))   # p-tic order (e.g., p = 6 for hexatic)
neigh_file = args.pop(1) # Neighbour list file
out_file = args.pop(1)   # Output file storing the order for each point

# Read the positions and neighbour lists
pos = []
neighbours = []
with open(neigh_file, 'r') as reader:
    for line in reader:
        if (line.startswith('#')): # Skip comments
            continue
        data = line.split()
        i = int(data[0])
        x = float(data[1])
        y = float(data[2])
        pos.append((x,y))
        nn = set()
        for j in range(3,len(data)):
            nn.add(int(data[j]))
        neighbours.append(nn)
pos = np.asarray(pos)
npoints = pos.shape[0]

# Compute and output the order parameter
def angle(x,y): # Get theta angle in the range [0,2*pi)
    return np.arctan2(-y,-x)+np.pi

with open(out_file, 'w') as writer:
    writer.write("#index x y Re(psi) Im(psi) |psi| Arg(psi)\n")    
    fmt = "{:d} {:g} {:g} {:g} {:g} {:g} {:g}\n"
    for i in range(npoints):
        psi = 0.0+0.0j
        nneigh = len(neighbours[i])
        for j in neighbours[i]:
            vec = pos[i]-pos[j]
            theta = angle(vec[0],vec[1])
            psi += np.exp(p*theta*1j)
        psi /= max(float(nneigh),1.0)
        writer.write(fmt.format(i, pos[i,0], pos[i,1], psi.real, psi.imag,
                                abs(psi), angle(psi.real,psi.imag)))
