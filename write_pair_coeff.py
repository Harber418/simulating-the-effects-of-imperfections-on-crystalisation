import numpy as np

# Prepare output lines for LJ/WCA
def Pair_coeffecients(epsilon=1, ntype=60,sizes=[],cut_factor=(2 ** (1 / 6)),filename = "pair.polydisperse"):
    pair_lines = []
    for i in range(ntype):
        for j in range(i, ntype):
            sigma = 0.5 * (sizes[i] + sizes[j])
            cutoff = cut_factor * sigma
            if filename == "pair.polydisperse.equilibration":
                pair_lines.append(f"pair_coeff {i+1} {j+1} {epsilon:.6f} {cutoff:.6f}")
            else:
                pair_lines.append(f"pair_coeff {i+1} {j+1} {epsilon:.6f} {sigma:.6f} {cutoff:.6f}")

    # Write all coefficients to file
    with open(filename, "w") as f:
        f.write("\n".join(pair_lines) + "\n")

    print(f" Pair coefficients written to {filename}")



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
    print(cut_factor)
    #produces equilibriation step 
    Pair_coeffecients(epsilon=100,sizes=size,filename="pair.polydisperse.equilibration")
    #Produces production step
    Pair_coeffecients(sizes=size)
    #equilibration for attractive step 
    Pair_coeffecients(sizes=size,cut_factor=1.5,filename = "pair.polydisperse.attractive.equilibration")
    #attractive step
    Pair_coeffecients(epsilon=5,sizes=size,cut_factor=1.5,filename = "pair.polydisperse.attractive")


if __name__ == "__main__":
    main()