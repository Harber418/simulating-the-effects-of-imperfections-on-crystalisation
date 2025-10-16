import numpy as np
import matplotlib.pyplot as plt

def generate_lammps(lx=35, ly=35, nbead=1000, seed=33,
                    mean=1.0, sd=0.1, filename="lammps_pd.data",
                    statsfile="type_stats.txt", make_plot=True):
    """
    Generate a 2D LAMMPS data file with random bead positions and
    Gaussian-distributed sizes binned into 60 discrete types.
    Also outputs type population statistics.
    """
    nbin = 60
    minsigma = 0.4
    maxsigma = 1.6
    binsize = (maxsigma - minsigma) / float(nbin)

    rng = np.random.default_rng(seed)

    # Random 2D positions
    #x = (rng.random(nbead) * 2 - 1) * lx / 2.0
    #y = (rng.random(nbead) * 2 - 1) * ly / 2.0
    #z = np.zeros(nbead)

    # Arrays
    bead_type = np.zeros(nbead, dtype=int)
    sigma = np.zeros(nbead)

    # Draw Gaussian sizes and assign to bins
    sizes_raw = rng.normal(loc=mean, scale=sd, size=nbead)
    for i in range(nbead):
        size = sizes_raw[i]
        if size <= minsigma:
            bead_type[i] = 1
        elif size >= maxsigma:
            bead_type[i] = nbin
        else:
            bead_type[i] = int((size - minsigma) / binsize) + 1
        sigma[i] = minsigma + (bead_type[i] - 0.5) * binsize

    # Compute distribution statistics
    mean_sigma = np.mean(sigma)
    std_sigma = np.std(sigma)
    PD = 100 * std_sigma / mean_sigma

    # Count how many beads of each type
    unique, counts = np.unique(bead_type, return_counts=True)
    type_counts = dict(zip(unique, counts))

    #fin the area fravtion 
    areaf = lx*ly/(1000*(np.pi*(mean/2)**2))
    print(f"area fraction initial is {areaf:.2f}")
    area_current = lx*ly/(np.sum(np.pi*(sigma/2)**2))

    lx = lx * np.sqrt(areaf/area_current)
    ly = ly * np.sqrt(areaf/area_current)
    print(f"updated box size to lx={lx:.2f}, ly={ly:.2f} to achieve area fraction {areaf:.3f}")
    #update the new size of the box 

    # Random 2D positions
    x = (rng.random(nbead) * 2 - 1) * lx / 2.0
    y = (rng.random(nbead) * 2 - 1) * ly / 2.0
    z = np.zeros(nbead)

    # --- Write LAMMPS data file ---
    with open(filename, "w") as f:
        f.write("LAMMPS 2D data file\n\n")
        f.write(f"{nbead} atoms\n")
        f.write(f"{nbin} atom types\n\n")
        f.write(f"{-lx/2:.6f} {lx/2:.6f} xlo xhi\n")
        f.write(f"{-ly/2:.6f} {ly/2:.6f} ylo yhi\n\n")

        # Masses section
        f.write("Masses\n\n")
        for i in range(nbin):
            mass_i = minsigma + (i + 0.5) * binsize
            f.write(f"{i+1} {mass_i:.6f}\n")
        f.write("\n")

        # Atoms
        f.write("Atoms # atomic\n\n")
        for i in range(nbead):
            atom_id = i + 1
            atom_type = bead_type[i]
            f.write(f"{atom_id} {atom_type} {x[i]:.6f} {y[i]:.6f} {z[i]:.6f}\n")
        f.write("\n")

        # Velocities (all zero)
        f.write("Velocities\n\n")
        for i in range(nbead):
            f.write(f"{i+1} 0.0 0.0 0.0\n")

    # --- Write statistics file ---
    with open(statsfile, "w") as f:
        f.write("Bead type statistics:\n")
        f.write("Type\tCount\n")
        for t in range(1, nbin + 1):
            count = type_counts.get(t, 0)
            f.write(f"{t}\t{count}\n")

        f.write("\n")
        f.write(f"Mean sigma = {mean_sigma:.6f}\n")
        f.write(f"Std sigma = {std_sigma:.6f}\n")
        f.write(f"Polydispersity (PD) = {PD:.2f}%\n")

    print(f"LAMMPS input written to '{filename}'")
    print(f"Type statistics written to '{statsfile}'")
    print(f"Mean sigma = {mean_sigma:.4f}, Std = {std_sigma:.4f}, PD = {PD:.2f}%")

    # Optional: plot histogram of types
    if make_plot:
        plt.figure(figsize=(8,4))
        plt.bar(unique, counts, width=0.8, color="skyblue", edgecolor="k")
        plt.xlabel("Bead Type (1–60)")
        plt.ylabel("Count")
        plt.title(f"Type Distribution (PD={PD:.2f}%)")
        plt.tight_layout()
        plt.show()

    return bead_type, sigma, type_counts


if __name__ == "__main__":
    generate_lammps(lx=35,ly=35, mean=1.0, sd=0.18, filename="lammps_input_pd", statsfile="type_stats.txt")