
"""
a tool to visualise the colloids in the box like vmd 
plots for colour reltates to size 
plot for colour relates to number of nearest neighbours
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

def read_better(file):
    coordinates=[]
    COLOUR=[]
    with open("dump.LJ", "r") as f:
        lines =f.readlines()
        for i in range(1,102):
            index1 = i*1009 +9 
            index2 = index1 + 1000
            set1 = []
            set2 = []
            atom_positions = lines[index1:index2]
            for atom in atom_positions:
                data =atom.split()
                c, x, y = float(data[1]),float(data[2]), float(data[3])
                set1.append([x,y])
                nbin = 60
                minsigma = 0.4
                maxsigma = 1.6
                binsize = (maxsigma - minsigma) / float(nbin)
                sigma = minsigma + (c - 0.5) * binsize
                set2.append(sigma)
            coordinates.append(np.array(set1, dtype=float))
            COLOUR.append(np.array(set2, dtype=float))

            
    return np.array(coordinates, dtype=object), COLOUR


def Number_near_neighbout_delinea(coordinates):
    tri = Delaunay(coordinates)
    triangles = tri.simplices
    thres =3.0

    def circumradius(try_pts):
        a = np.linalg.norm(try_pts[0]-try_pts[1])
        b = np.linalg.norm(try_pts[1]-try_pts[2])
        c = np.linalg.norm(try_pts[2]-try_pts[0])
        s = 0.5*(a+b+c)
        area = max(s*(s-a)*(s-b)*(s-c),0.0)**0.5
        if area == 0: return np.inf
        return (a*b*c)/(4.0*area)    
    # Perform Delaunay triangulation

    # Filter out triangles whose circumradius is too large
    triangles = tri.simplices
    tri_radii = np.asarray([circumradius(coordinates[t,:]) for t in triangles])
    med_radius = np.median(tri_radii)
    thres *= med_radius
    triangles = triangles[tri_radii < thres]

    # Create the neighbour lists - only consider the filtered triangles based on
    # the threshold on circumradius
    neighbours = [set() for i in range(len(coordinates))]
    for s in triangles:
        for i in s:
            for j in s:
                if (i != j):
                    neighbours[i].add(j)
    # Convert sets to sorted lists
    neighbours = [list(s) for s in neighbours]
    count = [len(s) for s in neighbours]
    return neighbours, count


def main():
    coordiantes, COLOUR = read_better("dump.LJ")
    snap= coordiantes[-1]
    Neighbors, count = Number_near_neighbout_delinea(snap)

    x,y = snap[:,0], snap[:,1]
    scatter_plot = plt.scatter(x , y,c=COLOUR[-1],marker='o',s=COLOUR[-1]*20)
    
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.title("scatter for last time step")
    cbar = plt.colorbar(scatter_plot)
    cbar.set_label("size of particles")

    plt.tight_layout()
    plt.show() 

def main2():
    coordiantes, COLOUR = read_better("dump.LJ")
    snap= coordiantes[-1]
    x,y = snap[:,0], snap[:,1]
    Neighbors, count = Number_near_neighbout_delinea(snap)
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.title("Neighbor count for crystal")
    sc = plt.scatter(x,y,c=count, marker='o', s=COLOUR[-1]*20)
    cbar = plt.colorbar(sc)
    cbar.set_label("number of nearest neighbours")
    plt.tight_layout()
    plt.show() 

if __name__ == "__main__":
    main()