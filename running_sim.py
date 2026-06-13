import subprocess

def run_lammps(filename="LJ2Dall.polydisperse.lam", log_file="log.lammps", lammps_exe="lmp", n_procs=1):
    cmd = ["mpirun", "-np", str(n_procs), lammps_exe, "-in", filename, "-log", log_file]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("LAMMPS run failed!")
        print(result.stderr)
        raise RuntimeError("LAMMPS simulation did not complete successfully")
    else:
        print("LAMMPS run completed successfully.")
    return result