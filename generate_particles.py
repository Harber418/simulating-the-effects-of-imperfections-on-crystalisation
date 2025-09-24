import numpy as np

def generate_lammps_input(lx=100.0, ly=100.0, ntype=1, nbead=1000, seed=71, filename="lammps_input"):
    """
    Generate a 2D LAMMPS data file with random bead positions.
    
    Parameters
    ----------
    lx, ly : float
        Box lengths in x and y directions
    ntype : int
        Number of atom types
    nbead : int
        Number of beads (atoms)
    seed : int
        Random number seed
    filename : str
        Output filename (default: lammps_input)
    """
    rng = np.random.default_rng(seed)

    # Random positions in [-lx/2, lx/2], [-ly/2, ly/2]
    x = (rng.random(nbead) * 2 - 1) * lx / 2.0
    y = (rng.random(nbead) * 2 - 1) * ly / 2.0
    z = np.zeros(nbead)

    with open(filename, "w") as f:
        f.write("LAMMPS 2D data file\n")
        f.write(f"{nbead} atoms\n\n")
        f.write(f"{ntype} atom types\n\n")
        f.write(f"{-lx/2:.6f} {lx/2:.6f} xlo xhi\n")
        f.write(f"{-ly/2:.6f} {ly/2:.6f} ylo yhi\n\n")

        f.write("Masses\n\n")
        for i in range(1, ntype + 1):
            f.write(f"{i} 1.0\n")

        f.write("\nAtoms # atomic\n\n")
        for i in range(nbead):
            atom_id = i + 1
            atom_type = 1
            f.write(f"{atom_id} {atom_type} {x[i]:.6f} {y[i]:.6f} {z[i]:.6f}\n")

        f.write("\nVelocities\n\n")
        for i in range(nbead):
            atom_id = i + 1
            f.write(f"{atom_id} 0.0 0.0 0.0\n")

    print(f"✅ LAMMPS input written to {filename}")


# Example usage:
if __name__ == "__main__":
    generate_lammps_input(lx=50.0, ly=50.0, ntype=1, nbead=100, seed=12345)
