import numpy as np

def write_lammps(attraction: bool = True, main_iterations: int = 1000000,
                  attraction_iterations: int = 250000, cut_off: float = 2.6,
                  timestep: float = 0.01):
    """
    bool : do you want attraction?
    int : iterations for main sequence
    int : iterations for attraction
    float : cut off
    float : timestep for attractive stage
    """
    filename = "LJ2Dall.polydisperse.lam"
    with open(filename, "w") as f:
        f.write("units lj\n")
        f.write("dimension 2 \n atom_style atomic \n")
        f.write("boundary p p p \n")
        f.write("neighbor 1.9 bin\n")
        f.write("neigh_modify every 1 delay 1 check yes \n")
        f.write("read_data lammps_input_pd\n")
        f.write("reset_timestep 0\n")

        # equilibration starts here
        f.write("pair_style soft 1.12246152962189\n")
        f.write("include pair.polydisperse.equilibration\n")
        #fix langevin Tstart Tstop Tdamp seed
        f.write("fix 1 all nve\nfix 2 all langevin 1.0 1.0 1.0 31\nfix 3 all enforce2d\n")
        f.write("thermo 1000\n")
        f.write("thermo_style custom step temp epair press pxx pyy vol\n")
        f.write("timestep 0.01\nrun 250000\n")
        f.write("write_dump all atom equil_positions.dump\n")

        # production stage
        f.write("pair_style lj/cut 1.12246152962189\n")
        f.write("pair_modify shift yes\n")
        f.write("include pair.polydisperse\n")
        f.write("thermo 10000\n")
        f.write("thermo_style custom step temp epair press pxx pyy vol\n")
        f.write("timestep 0.01\n")
        f.write(f"run {main_iterations}\n")
        f.write("write_dump all atom production_positions.dump\n")

        if attraction:
            # equilibration for attractive stage
            f.write(f"pair_style lj/cut {cut_off}\n")
            f.write("pair_modify shift no\n")
            f.write("include pair.polydisperse.attractive.equilibration\n")
            f.write(f"timestep {timestep}\nrun 100000\n")
            f.write("write_dump all atom attractive_equil_positions.dump\n")

            # attractive stage
            f.write(f"pair_style lj/cut {cut_off}\npair_modify shift no\ninclude pair.polydisperse.attractive\n")
            f.write(f"thermo 10000\nthermo_style custom step temp epair press pxx pyy vol\ntimestep {timestep}\n")
            f.write(f"run {attraction_iterations}\n")
            f.write("write_dump all atom attractive_positions.dump\n")