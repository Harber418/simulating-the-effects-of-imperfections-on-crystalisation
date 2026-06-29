"""
This file is used to run the whole simulation pipeline 
from given inputs is 
genetates particle postions for desired polydispersity 
writes the pair correlation focuntion for the different particle interactions 
it writes the lammps file with desiered attraction parameters 
runs the lammps simulation 
saves the required data as a text file 
produces plots of the psi 6 for each step in the simulation.
"""
#https://iifx.dev/en/articles/457242272/bridging-the-gap-how-to-make-your-c-matrix-multiplications-as-fast-as-fortran
import numpy as np 
import sys
import argparse
from Generate_lammps_area_corrected import generate_lammps
from write_pair_coeff import Pair_coeffecients
from write_lammps_input_file import write_lammps
from running_sim import run_lammps
from visualise import organise_and_plot

def main():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-L', '--size', type=float, default=30, help='System size (default: 30)')
    parser.add_argument('-P', '--polydispersity', type=float, default=9.6, help='polydispersity (default: 9%)')
    parser.add_argument('-E', '--epsilon', type=float, default=1, help='attractive constant (default: 1)')
    parser.add_argument('-c', '--cut_off', type=float, default=1.05, help="size of cutoff for attractive tail (default: '1.5')")
    parser.add_argument('-a', '--attractive', type=bool, default=True, help="do you want to add the attractive stage to the simulation (default: True)")
    parser.add_argument('-i', '--iterations', type=int, default=1500000, help="how long should the production stage be (default: 1000000)")
    parser.add_argument('-I', '--attractive_itterations', type=int, default=100000, help="how long should the attractive production stage be (default: 250000)")
    parser.add_argument('-t', '--timestep', type=float, default=0.01, help="timestep for attractive term, used if forces are too large leading to errors (default: 0.01)")
    parser.add_argument('-s', '--scaled' ,type=bool,default=True, help = "diffenert sized particles attract differently ")
    parser.add_argument('-v', '--volume' ,type=float,default=32, help = "length of box to shrink to")
    args = parser.parse_args() 


    #===========================================================================================
    #generate initial positions for particles and give sizes 
    generate_lammps(lx=args.size,ly=args.size, mean=1.0, sd=args.polydispersity*0.01, filename="lammps_input_pd", statsfile="type_stats.txt")

    #===========================================================================================================
    #write pair coeffecients 
    ntype = 60
    minsigma = 0.4
    maxsigma = 1.6

    # Compute size array
    binsize = (maxsigma - minsigma) / ntype
    size = np.array([minsigma + (i + 0.5) * binsize for i in range(ntype)])

    # WCA cutoff factor
    cut_factor = 2 ** (1 / 6) #1.12246
    print(cut_factor)
    #produces equilibriation step 
    Pair_coeffecients(epsilon=100,sizes=size,filename="pair.polydisperse.equilibration")
    #Produces production step
    Pair_coeffecients(epsilon=10, sizes=size)
    #equilibration for attractive step 
    Pair_coeffecients(epsilon=args.epsilon,sizes=size,cut_factor=args.cut_off,filename = "pair.polydisperse.attractive.equilibration",scaled=args.scaled)
    #attractive step
    Pair_coeffecients(epsilon=args.epsilon,sizes=size,cut_factor=args.cut_off,filename = "pair.polydisperse.attractive",scaled=args.scaled)
    
    
    #==============================================================================================================================
    #write the lammps file 
    write_lammps(attraction=args.attractive,main_iterations=args.iterations,attraction_iterations=args.attractive_itterations,cut_off=args.cut_off*1.5,timestep=args.timestep,volume_reduction=args.volume)
    
    #================================================================================================================================
    #run the lammps simulation 
    print("running lammps simulation please wait")
    print("this may take up to 10 minutes depending on choosen parameters")
    run_lammps(filename="LJ2Dall.polydisperse.lam", lammps_exe="lmp")

    #=================================================================================
    #now we want to save the images for each stage of production to add validity to the process
    if args.attractive:
        organise_and_plot(output_dir=f"results_PD{args.polydispersity}_size{args.size}_attractive_epsilon{args.epsilon}_cutoff{args.cut_off}")
    else:
        organise_and_plot(output_dir=f"results_PD{args.polydispersity}_size{args.size}")


if __name__ == "__main__":
    main()