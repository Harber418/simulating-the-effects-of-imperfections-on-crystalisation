
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

# Prepare output lines for LJ/WCA
pair_lines = []
for i in range(ntype):
    for j in range(i, ntype):
        epsilon = 1.0
        sigma = 0.5 * (sizes[i] + sizes[j])
        cutoff = cut_factor * sigma
        pair_lines.append(f"pair_coeff {i+1} {j+1} {epsilon:.6f} {sigma:.6f} {cutoff:.6f}")

# Write all coefficients to file
with open("pair.polydisperse", "w") as f:
    f.write("\n".join(pair_lines) + "\n")

print(" Pair coefficients written to 'pair.polydisperse'")

# Prepare output lines for soft
pair_lines = []
for i in range(ntype):
    for j in range(i, ntype):
        epsilon = 100.0
        sigma = 0.5 * (sizes[i] + sizes[j])
        cutoff = cut_factor * sigma
        pair_lines.append(f"pair_coeff {i+1} {j+1} {epsilon:.6f} {cutoff:.6f}")

# Write all coefficients to file
with open("pair.polydisperse.equilibration", "w") as f:
    f.write("\n".join(pair_lines) + "\n")

print(" Pair coefficients written to 'pair.polydisperse.equilibration'")