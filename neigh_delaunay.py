#!/usr/bin/env python3
# neigh_delaunay.py
# A script to find the nearest neighbours using Delaunay triangulation.

import sys
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay

args = sys.argv
if (len(args) != 5):
    print("Usage: neigh_delaunay.py xcol ycol pos_file out_file")
    sys.exit(1)

xcol = int(args.pop(1)) # Column index of the x position of points
ycol = int(args.pop(1)) # Column index of the y position of points
pos_file = args.pop(1)  # Position file
out_file = args.pop(1)  # Output file storing the neighbour list for each point

# Threshold (in units of the median circumradius) for determining if a triangle
# should be kept when doing Delaunay triangulation. This tries to eliminate
# (most of) the extraneous neighbours at the boundaries
thres = 3.0

# For computing the radius of the circumcircle of each triangle
def circumradius(tri_pts):
    a = np.linalg.norm(tri_pts[0]-tri_pts[1])
    b = np.linalg.norm(tri_pts[1]-tri_pts[2])
    c = np.linalg.norm(tri_pts[2]-tri_pts[0])
    s = (a+b+c)*0.5
    area = max(s*(s-a)*(s-b)*(s-c),0.0)**0.5
    if (area == 0):
        return np.inf
    return (a*b*c)/(4.0*area)

# Read the position data
pos = pd.read_csv(pos_file, delim_whitespace=True,usecols=[xcol,ycol]).to_numpy()
npoints = pos.shape[0]

# Perform Delaunay triangulation
tri = Delaunay(pos)

# Filter out triangles whose circumradius is too large
triangles = tri.simplices
tri_radii = np.asarray([circumradius(pos[t,:]) for t in triangles])
med_radius = np.median(tri_radii)
thres *= med_radius
triangles = triangles[tri_radii < thres]

# Create the neighbour lists - only consider the filtered triangles based on
# the threshold on circumradius
neighbours = [set() for i in range(len(pos))]
for s in triangles:
    for i in s:
        for j in s:
            if (i != j):
                neighbours[i].add(j)

# Output the neighbour list for each point (only those within the box but not
# the ones in the buffered region)
with open(out_file, 'w') as writer:
    writer.write("#index x y neighbours\n")
    for i in range(npoints):
        # For convenience, output the index and position for each point
        writer.write("{:d} {:g} {:g} ".format(i, pos[i][0], pos[i][1]))
        for n in sorted(neighbours[i]):
            writer.write(f"{n} ")
        writer.write("\n")
