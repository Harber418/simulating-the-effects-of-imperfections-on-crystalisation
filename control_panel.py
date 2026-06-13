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

from Generate_lammps_area_corrected import generate_lammps
from write_pair_coeff import Pair_coeffecients
from write_lammps_input_file import write_lammps
