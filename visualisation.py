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
                b = 0.1
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

def angles_for_NN(coordinates, neighbors):
    angles =[]
    for i in range(len(coordinates)):
        #for each atom find the 5 angles that make up the hexagonal strucutre and add it to a list  
        #atom i
        xi, yi = coordinates[i]
        # near neiborours for atom i (already sorted)
        Nn = np.array(coordinates[neighbors[i]], dtype=float)
        #compure vecotrs from atom i to near neibours 
        #and normalise 
        vectors = Nn - np.array([xi,yi], dtype=float)
        norms = np.linalg.norm(vectors,axis=1, keepdims=True)
        norms[norms == 0] = 1.0 
        vectors = vectors /norms

        #now compute all 6 anlges around one atom  
        angles_for_atomi =[]
        for A in vectors:
            theta = np.arctan2(A[1], A[0]) 
            angles_for_atomi.append(theta)
        # append atom i's 6 angle to a master list for all atoms 
        angles.append(angles_for_atomi)
    return angles


def psi_six_local_order(angles):
    "takes angles and nearest neighbour and returns psi 6 value"
    psi = []
    for angle_list in angles:
        psi_i = 0 
        for x in angle_list: 
            exponential = np.exp(6j*x)
            psi_i += exponential
        psi.append(psi_i/len(angle_list))

    return np.array(psi)

def psi_six_global(local,N):
    "golbal version "
    psi = np.abs(np.mean(local))
    return psi

def psi_plot():
    coordiantes, COLOUR = read_better("dump.LJ")
    snap= coordiantes[-1]
    nn, cc = Number_near_neighbout_delinea(snap)
    anlge_input = angles_for_NN(snap, nn)
    N = len(snap)
    psi= psi_six_local_order(anlge_input)
    global_psi = psi_six_global(psi,N)
    x,y = snap[:,0], snap[:,1]
    mean_sigma = np.mean(COLOUR[-1])
    std_sigma = np.std(COLOUR[-1])
    PD = 100 * std_sigma / mean_sigma
    scatter_plot = plt.scatter(x , y,c=psi,marker='o',cmap='jet',s=COLOUR[-1]*20, vmin=0.45, vmax=0.93)
    plt.xlabel("position x ")
    plt.ylabel("position y ")
    plt.title(f"Colloid crystal Psi6 global = {global_psi:.3f} (polydispersity {PD:.2f}%)")
    cbar = plt.colorbar(scatter_plot)
    cbar.set_label("size of particles")

    plt.tight_layout()
    plt.show() 




def size():
    coordiantes, COLOUR = read_better("dump.LJ")
    snap= coordiantes[-1]


    x,y = snap[:,0], snap[:,1]
    mean_sigma = np.mean(COLOUR[-1])
    std_sigma = np.std(COLOUR[-1])
    PD = 100 * std_sigma / mean_sigma
    scatter_plot = plt.scatter(x , y,c=COLOUR[-1],marker='o',cmap='jet',s=COLOUR[-1]*20, vmin=0.54, vmax=1.4)
    
    plt.xlabel("position x ")
    plt.ylabel("position y ")
    plt.title(f"Colloid crystal (polydispersity {PD:.2f}%)")
    cbar = plt.colorbar(scatter_plot)
    cbar.set_label("size of particles")

    plt.tight_layout()
    plt.show() 
    
def main():
    psi_plot()
    size()
    
    coordiantes, COLOUR = read_better("dump.LJ")
    snap= coordiantes[-1]
    x,y = snap[:,0], snap[:,1]
    Neighbors, count = Number_near_neighbout_delinea(snap)
    b=0.05
    mask = (x >= b) & (x <= 1 - b) & (y >= b) & (y <= 1 - b)

    # Apply mask to get "interior" particles
    X = x[mask]
    Y = y[mask]
    counts_inside = np.array(count)[mask]
    color_inside = COLOUR[-1][mask]
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.title("Neighbor count for crystal")
    sc = plt.scatter(X,Y,c=counts_inside, marker='o', s=color_inside*20,vmin=2, vmax=8)
    cbar = plt.colorbar(sc)
    cbar.set_label("number of nearest neighbours")
    plt.tight_layout()
    plt.show() 
    
if __name__ == "__main__":
    main()