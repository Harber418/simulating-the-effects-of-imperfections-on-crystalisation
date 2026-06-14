"""
a tool to visualise the colloids in the box like vmd 
plots for colour reltates to size 
plot for colour relates to number of nearest neighbours
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import cm, colors
from scipy.spatial import Delaunay
import os
import shutil
import scienceplots
plt.style.use('science') # more scientific style for matplotlib
plt.rcParams['text.usetex'] = False # this avoids an annoying latex installation


def read_last(file):
    
    nbin = 60
    minsigma = 0.4
    maxsigma = 1.6
    binsize = (maxsigma - minsigma) / float(nbin)

    with open(file, "r") as f:
        lines = f.readlines()

    # Find the LAST occurrence of "ITEM: ATOMS"
    last_idx = None
    box_idx = None
    for i in range(len(lines) - 1, -1, -1):
        if last_idx is None and lines[i].strip() == "ITEM: ATOMS id type xs ys zs":
            last_idx = i + 1
        if box_idx is None and lines[i].startswith("ITEM: BOX BOUNDS"):
            box_idx = i + 1
        if last_idx is not None and box_idx is not None:
            break

    if last_idx is None or box_idx is None:
        raise ValueError("No atom/box data found in dump file")

    xlength, xhight = map(float, lines[box_idx].split()[:2])
    #line looks like this -1.4916282000000001e+01 1.4916282000000001e+01
    Lx = np.abs(xhight - xlength)  # box is square, so this also equals Ly


    set1 = []
    set2 = []
    i = last_idx
    while i < len(lines) and not lines[i].startswith("ITEM:"):
        data = lines[i].split()
        c = float(data[1])
        x = float(data[2])
        y = float(data[3])
        sigma = minsigma + (c - 0.5) * binsize
        set1.append([x, y])
        set2.append(sigma/ Lx) # we normalise sigma so it fits the coordiantes scale between 0 and 1 from lammps
        i += 1

    return np.array(set1, dtype=float), np.array(set2, dtype=float)


def save_data(file):
    snap, COLOUR = read_last("dump.LJ")
    #we want to save the last data so that we can read it later
    #as each dump.LJ is 5 MB we want to reduce this for larger data collection. 
    with open(file, "a") as f:
        f.write("sigma, x , y\n")
        for step, (M, E) in enumerate(zip(snap, COLOUR)):
            f.write(f"{E},{M[0]},{M[1]}\n")



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
    #print(f"the threashhold for delaunay is = {thres}")
    triangles = triangles[tri_radii < thres]
    #print(triangles)

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
    psi = np.mean(np.abs(local))
    return psi

def psi_plot(filename="dump.LJ", save_path=None,Phase = True):
    snap, COLOUR = read_last(filename)
    nn, cc = Number_near_neighbout_delinea(snap)
    angle_input = angles_for_NN(snap, nn)
    N = len(snap)

    psi = psi_six_local_order(angle_input)
    #global_psi = psi_six_global(psi, N)
    x, y = snap[:, 0], snap[:, 1]

    mean_sigma = np.mean(COLOUR)
    std_sigma = np.std(COLOUR)
    PD = 100 * std_sigma / mean_sigma

    b = 0.05
    mask = (x >= b) & (x <= 1 - b) & (y >= b) & (y <= 1 - b)
    X, Y = x[mask], y[mask]
    psi = psi[mask]
    color_inside = COLOUR[mask]
    global_psi = np.mean(np.abs(psi))#global claculated based on the internal sample
    fig, ax = plt.subplots()

    cmap = cm.get_cmap('Wistia')
    norm = colors.Normalize(vmin=0.0, vmax=1.0)

    for xi, yi, sigma_i, psi_i in zip(X, Y, color_inside, psi):
        color = cmap(norm(np.abs(psi_i)))
        circ = Circle((xi, yi), radius=sigma_i/2, facecolor=color, edgecolor='k',linewidth=0.1)
        ax.add_patch(circ)
    ax.set_aspect("equal")
    # add colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([]) 
    cbar = fig.colorbar(sm,ax=ax)
    cbar.set_label(r"$|\psi_6|$")
    #fig.colorbar(sm, ax=ax, label='Local |psi6|')
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.text(0.02,0.98,rf"$Global \phi={global_psi:.3f}$"+"\n"+rf"$Polydispersity={PD:.2f}\%$",transform=ax.transAxes,va="top")
    ax.set_title(rf"$Global " +" "+rf"\phi:{global_psi:.3f}$"+ " "+rf"$Polydispersity:{PD:.2f}\%$")
    
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        #colloidal_distance(snap)
    else:
        plt.show()

    #
    #add a phase diagram 
    if Phase:
        fig, ax = plt.subplots()
        phase_norm = colors.Normalize(vmin=-np.pi,vmax=np.pi)
        cmap = cm.get_cmap('turbo')
        for xi, yi, sigma_i, psi_i in zip(X, Y, color_inside, psi):
            phase = np.angle(psi_i)/6
            color = plt.cm.hsv(phase_norm(phase))
            circ = Circle((xi, yi), radius=sigma_i/2, facecolor=color, edgecolor='k',linewidth=0.1)
            ax.add_patch(circ)
        ax.set_aspect("equal")
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])

        sm = cm.ScalarMappable(
            cmap=plt.cm.hsv,
            norm=phase_norm
        )
        sm.set_array([])

        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label(r'arg($\psi_6$)')

        ax.set_title(
            rf"$Global\ \phi={global_psi:.3f}$"
            + " "
            + rf"$PD={PD:.2f}\%$"
        )

        plt.tight_layout()
        root, ext = os.path.splitext(save_path)

        phase_path = root + "_phase" + ext

        fig.savefig(
            phase_path,
            dpi=300,
            bbox_inches="tight"
        )

        print(f"Saved phase plot: {phase_path}")

        plt.close(fig)

def psi_plot_old(filename="dump.LJ", save_path=None):
    snap, COLOUR = read_last(filename)
    nn, cc = Number_near_neighbout_delinea(snap)
    angle_input = angles_for_NN(snap, nn)
    N = len(snap)

    psi = psi_six_local_order(angle_input)
    global_psi = psi_six_global(psi, N)
    x, y = snap[:, 0], snap[:, 1]

    mean_sigma = np.mean(COLOUR)
    std_sigma = np.std(COLOUR)
    PD = 100 * std_sigma / mean_sigma

    b = 0.05
    mask = (x >= b) & (x <= 1 - b) & (y >= b) & (y <= 1 - b)
    X, Y = x[mask], y[mask]
    psi = psi[mask]
    color_inside = COLOUR[mask]

    fig, ax = plt.subplots()
    #scatter_plot = ax.scatter(X, Y, c=np.abs(psi), marker='o', cmap='Wistia',
    #                           s=color_inside*20, vmin=0.3, vmax=0.93)
    ax.set_xlabel("position x")
    ax.set_ylabel("position y")
    ax.set_title(f"Psi6 global = {global_psi:.3f} (polydispersity {PD:.2f}%)")
    #cbar = fig.colorbar(scatter_plot)
    #cbar.set_label("Local psi 6")
    plt.tight_layout()
    #new plotting 
    #fig, ax = plt.subplots()

    for x,y,sigma in zip(X,Y,COLOUR):
        circle = Circle(
            (x,y),
            radius=sigma/2,
            fill=True,
            c=np.abs(psi),
            cmap='Wistia'
        )
        ax.add_patch(circle)

    ax.set_aspect("equal")

    #

    if save_path:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        #colloidal_distance(snap)
    else:
        plt.show()
    


def organise_and_plot(dump_files=None, output_dir="results"):
    """
    Moves the given dump files into output_dir and saves a psi6 plot for each.
    """
    if dump_files is None:
        dump_files = [
            "equil_positions.dump",
            "production_positions.dump",
            "attractive_equil_positions.dump",
            "attractive_positions.dump",
        ]

    os.makedirs(output_dir, exist_ok=True)

    for dump_file in dump_files:
        if not os.path.exists(dump_file):
            print(f"Note: {dump_file} was skiped.")
            continue

        # Generate the psi6 plot BEFORE moving the file
        stage_name = os.path.splitext(dump_file)[0]
        plot_path = os.path.join(output_dir, f"{stage_name}_psi6.png")
        psi_plot(filename=dump_file, save_path=plot_path)
        print(f"Saved plot: {plot_path}")

        # Move the dump file into the results folder
        dest = os.path.join(output_dir, dump_file)
        shutil.move(dump_file, dest)
        print(f"Moved {dump_file} to {dest}")
        #move phase 
        root, ext = os.path.splitext(plot_path)
        phase_path = root + "_phase" + ext
        phase_dest = os.path.join(output_dir, os.path.basename(phase_path))
        shutil.move(phase_path, phase_dest)



def colloidal_distance(coordinates):

    neighbours, count = Number_near_neighbout_delinea(coordinates)

    distances = []

    for i, neighs in enumerate(neighbours):
        for j in neighs:
            if j > i:  # avoid double counting
                d = np.linalg.norm(coordinates[i] - coordinates[j])
                distances.append(d)

    distances = np.array(distances)

    plt.hist(distances, bins=100)
    plt.xlabel("Neighbour distance")
    plt.ylabel("Count")
    plt.show()

def size():
    #coordiantes, COLOUR = read_better("dump.LJ")
    #snap= coordiantes[-1]
    snap, COLOUR = read_last("dump.LJ")

    x,y = snap[:,0], snap[:,1]
    mean_sigma = np.mean(COLOUR)
    std_sigma = np.std(COLOUR)
    PD = 100 * std_sigma / mean_sigma
    scatter_plot = plt.scatter(x , y,c=COLOUR,marker='o',cmap='jet',s=COLOUR[-1]*20, vmin=0.54, vmax=1.4)
    
    plt.xlabel("position x ")
    plt.ylabel("position y ")
    plt.title(f"Colloid crystal (polydispersity {PD:.2f}%)")
    cbar = plt.colorbar(scatter_plot)
    cbar.set_label("size of particles")

    plt.tight_layout()
    plt.show() 
    
def main():
    filename= "attraction.txt" 
    save_data(filename) #here we generate the data from a lammps file 
    #so if you want to look at old data.
    psi_plot()
    size()
    
    #coordiantes, COLOUR = read_better("dump.LJ")
    #snap= coordiantes[-1]
    snap, COLOUR = read_last("dump.LJ")
    x,y = snap[:,0], snap[:,1]
    Neighbors, count = Number_near_neighbout_delinea(snap)
    b=0.05
    mask = (x >= b) & (x <= 1 - b) & (y >= b) & (y <= 1 - b)

    # Apply mask to get "interior" particles
    X = x[mask]
    Y = y[mask]
    counts_inside = np.array(count)[mask]
    color_inside = COLOUR[mask]
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.title("Neighbor count for crystal")
    sc = plt.scatter(X,Y,c=counts_inside, marker='o', s=color_inside*20,vmin=3, vmax=8)
    cbar = plt.colorbar(sc)
    cbar.set_label("number of nearest neighbours")
    plt.tight_layout()
    plt.show() 
    
if __name__ == "__main__":
    main()
