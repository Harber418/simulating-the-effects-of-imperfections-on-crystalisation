import numpy as np

# Parameters
ntype = 60
minsigma = 0.4
maxsigma = 1.6

# Compute size array
binsize = (maxsigma - minsigma) / ntype
sizes = np.array([minsigma + (i + 0.5) * binsize for i in range(ntype)])

# WCA cutoff factor
cut_factor = 2 ** (1 / 6)
print(cut_factor)
# Prepare output lines for LJ/WCA
def Pair_coeffecients(epsilon=1, ntype=60,sizes=[],cut_factor=(2 ** (1 / 6)),filename = "pair.polydisperse"):
    pair_lines = []
    for i in range(ntype):
        for j in range(i, ntype):
            sigma = 0.5 * (sizes[i] + sizes[j])
            cutoff = cut_factor * sigma
            pair_lines.append(f"pair_coeff {i+1} {j+1} {epsilon:.6f} {sigma:.6f} {cutoff:.6f}")

    # Write all coefficients to file
    with open(filename, "w") as f:
        f.write("\n".join(pair_lines) + "\n")

    print(f" Pair coefficients written to {filename}")

# Prepare output lines for soft
pair_lines = []
for i in range(ntype):
    for j in range(i, ntype):
        epsilon = 100.0 #this is overwriten in the lammps simulation as we set pair_style soft 1.12246152962189
        sigma = 0.5 * (sizes[i] + sizes[j])
        cutoff = cut_factor * sigma
        pair_lines.append(f"pair_coeff {i+1} {j+1} {epsilon:.6f} {cutoff:.6f}")

# Write all coefficients to file
with open("pair.polydisperse.equilibration", "w") as f:
    f.write("\n".join(pair_lines) + "\n")

print(" Pair coefficients written to 'pair.polydisperse.equilibration'")

#Writing the coefficients for the attractive stage. 
#we need to vary epsilon but reasonable values should be 1 and 5 

epsilon_attractive = 10.0   # EXPERIMENTAL PARAMETER
attractive_cutoff_factor = 1.3   # extend beyond the WCA point to include the well

pair_lines = []
for i in range(ntype):
    for j in range(i, ntype):
        sigma = 0.5 * (sizes[i] + sizes[j])
        cutoff = attractive_cutoff_factor * sigma
        pair_lines.append(
            f"pair_coeff {i+1} {j+1} {epsilon_attractive:.6f} {sigma:.6f} {cutoff:.6f}"
        )

with open("pair.polydisperse.attractive", "w") as f:
    f.write("\n".join(pair_lines) + "\n")

print(f"Attractive pair coefficients written (epsilon={epsilon_attractive})")


def main():
    # Parameters
    ntype = 60
    minsigma = 0.4
    maxsigma = 1.6

    # Compute size array
    binsize = (maxsigma - minsigma) / ntype
    size = np.array([minsigma + (i + 0.5) * binsize for i in range(ntype)])

    # WCA cutoff factor
    cut_factor = 2 ** (1 / 6)
    #produces equilibriation step 
    Pair_coeffecients(epsilon=100,sizes=size,filename="pair.polydisperse.equilibration")
    #Produces production step
    Pair_coeffecients(sizes=size)
    #equilibration for attractive step 
    Pair_coeffecients(sizes=size,cut_factor=1.3,filename = "pair.polydisperse.attractive.equilibration")
    #attractive step
    Pair_coeffecients(epsilon=10,sizes=size,cut_factor=1.3,filename = "pair.polydisperse.attractive")
