
"""
Code that pulls lammps input sizes from a gaussian 
"""
import numpy as np 

def generate_lammps(lx=30,ly=30,ntype=100,nbead=1000, seed = 33,mean=1,sd=0.1, filename="lamps_pd"):
    """
    Generate a 2D LAMMPS data file with random bead positions and masses 
    pulled from a gaussian distribution
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
        masses = []
        for i in range(1, nbead + 1):
            # each atom gets a unique mass pulled form gaussian
            m = rng.normal(loc=mean, scale=sd)  # Mean 1.0, stddev 0.1
            # then each mass is binned into nytype bins
            masses.append(m)
        #make nypte evenly spaced bins between max and min values 
        bins = np.linspace(min(masses), max(masses), ntype+1)
        inds = np.digitize(masses, bins)
         #for each bin write the masses to a list 
        bin_full =[]
        for i in range(1, ntype + 1):
            bin_full.append(np.sum(inds==i))

        #find the center of each bin
        # low edge + high edge and dive by 2 for mid point 
        center = 0.5 * (bins[:-1] + bins[1:])
        #write masses as atom types to file 
        for i in range(1, ntype +1):
            f.write(f"{i} {center[i-1]}\n")
        

        f.write("\nAtoms # atomic\n\n")
        for i in range(nbead):
            atom_id = i + 1
            atom_type = inds[1]
            f.write(f"{atom_id} {atom_type} {x[i]:.6f} {y[i]:.6f} {z[i]}\n")

        f.write("\nVelocities\n\n")
        for i in range(nbead):
            atom_id = i + 1
            f.write(f"{atom_id} 0.0 0.0 0.0\n")

    print(f"LAMMPS input written to {filename} with 100 bins")
    sigma = np.std(masses)
    sample_mean = np.mean(masses)
    
    PD = 100 * sigma/sample_mean
    print(f"The polydispericty is {PD:.6f}%with a sample sigma of {sigma:.6f} and a sample mean of {sample_mean:.4f}")
    return 


if __name__ == "__main__":
    generate_lammps(ntype=100, mean=1.0, sd=0.1, filename="lammps_pd")