import os
import sys
import src.vars

dirname = os.path.dirname
path = os.path.dirname(__file__)

def get_swift_flavour(hydro, kernel):
    
    assert hydro in ("nSPH", "GDF"), f"hydro must be nSPH, or GDF"
    assert kernel in ("s3", "wc6"), f"kernel must be s3, or wc6"
    
    if hydro == "nSPH" and kernel == "s3":
        swift = "/cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/swift_nSPH_s3"
    elif hydro == "nSPH" and kernel == "wc6":
        swift = "/cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/swift_nSPH_wc6"
    elif hydro == "GDF" and kernel == "s3":
        swift = "/cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/swift_GDF_s3"
    elif hydro == "GDF" and kernel == "wc6":
        swift = "/cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/swift_GDF_wc6"
    else:
        raise NotImplementedError(f"Combination of hydro scheme and kernel not compiled yet.")
        
    return swift

def get_dir_name(N, m, **kwargs):
    
    L = kwargs.get("L", 0)
    kernel = kwargs.get("kernel", "wc6")
    hydro = kwargs.get("hydro", "GDF")
    
    assert N in (1e5, 1e6, 1e7), f"N must be 1e5, 1e6, or 1e7"
    assert m > 0 and m < 1e-8*src.vars.M_earth, f"m must be in M_earth units"
    assert hydro in ("nSPH", "GDF"), f"hydro must be nSPH, or GDF"
    assert kernel in ("s3", "wc6"), f"kernel must be s3, or wc6"
    
    name = "m" + str(m).replace(".", "") + "_L" + str(L).replace(".", "") + "_hydro" + hydro + "_kernel" + kernel
    return name

def get_dir_settling(N, m, **kwargs):
    
    assert N in (1e5, 1e6, 1e7), f"N must be 1e5, 1e6, or 1e7"
    assert m < 1e-8*src.vars.M_earth, f"m must be in M_earth units"
    
    N_str = "N{:.1e}".format(N).replace("+", "").replace(".0", "").replace("e0", "e")
    name = get_dir_name(N, m, **kwargs)
    folder_name = os.path.join(dirname(path), "data/sims", N_str, name)
    
    return folder_name

def make_dir_settling(N, m, **kwargs):
    
    folder_name = get_dir_settling(N, m, **kwargs)
    
    if not os.path.isdir(folder_name):
        command = "mkdir " + folder_name
        os.system(command)
        
        command = "mkdir " + os.path.join(folder_name, "snapshots")
        os.system(command)
        command = "mkdir " + os.path.join(folder_name, "restart")
        os.system(command)
        command = "mkdir " + os.path.join(folder_name, "images")
        os.system(command)
        
    else:
        print(f"Folder {folder_name} already exists.")
    
    return None

def make_yml_settling(N, m, **kwargs):
    
    delta_time = kwargs.get("delta_time", 2000)
    eta = kwargs.get("eta", 1.2348)
    h_max = kwargs.get("h_max", 0.1)
    mat = kwargs.get("mat", "ANEOS")
    
    assert mat in ("ANEOS", "Til"), f"mat must be in ANEOS or Til"
    
    folder_name = get_dir_settling(N, m, **kwargs)
    fn = os.path.join(folder_name, "planet.yml")
    
    f = open(fn,"w+")

    f.write("# Define the system of units to use internally.\n")
    f.write("InternalUnitSystem:\n")
    f.write("    UnitMass_in_cgs:        5.9724e27   # Grams\n")
    f.write("    UnitLength_in_cgs:      6.371e8     # Centimeters\n")
    f.write("    UnitVelocity_in_cgs:    6.371e8     # Centimeters per second\n")
    f.write("    UnitCurrent_in_cgs:     1           # Amperes\n")
    f.write("    UnitTemp_in_cgs:        1           # Kelvin\n\n")

    f.write("# Parameters related to the initial conditions\n")
    f.write("InitialConditions:\n")
    f.write("    file_name:  ./planet.hdf5   # The initial conditions file to read\n")
    f.write("    periodic:   0               # Are we running with periodic ICs?\n\n")

    f.write("# Parameters governing the time integration\n")
    f.write("TimeIntegration:\n")
    f.write("    time_begin:     0                   # The starting time of the simulation (in internal units).\n")
    f.write("    time_end:       20000              # The end time of the simulation (in internal units).\n")
    f.write("    dt_min:         0.0001              # The minimal time-step size of the simulation (in internal units).\n")
    f.write("    dt_max:         1000                # The maximal time-step size of the simulation (in internal units).\n\n")

    f.write("# Parameters governing the snapshots\n")
    f.write("Snapshots:\n")
    f.write("    basename:       ./snapshots/planet  # Common part of the name of output files\n")
    f.write("#    output_list_on: 0\n")
    f.write("#    output_list:    output_list.txt\n")
    f.write("    time_first:     0                   # Time of the first output (in internal units)\n")
    f.write("    delta_time:     " + "{:d}".format(delta_time) + "                # Time difference between consecutive outputs (in internal units)\n")
    f.write("    int_time_label_on:  1               # Enable to label the snapshots using the time rounded to an integer (in internal units)\n\n")

    f.write("# Parameters governing the conserved quantities statistics\n")
    f.write("Statistics:\n")
    f.write("    time_first:     0                   # Time of the first output (in internal units)\n")
    f.write("    delta_time:     1000                # Time between statistics output\n\n")

    f.write("# Parameters controlling restarts\n")
    f.write("Restarts:\n")
    f.write("    enable:         1                   # Whether to enable dumping restarts at fixed intervals.\n")
    f.write("    save:           1                   # Whether to save copies of the previous set of restart files (named .prev)\n")
    f.write("    onexit:         1                   # Whether to dump restarts on exit (*needs enable*)\n")
    f.write("    subdir:         ./restart           # Name of subdirectory for restart files.\n")
    f.write("    basename:       planet              # Prefix used in naming restart files.\n")
    f.write("    delta_hours:    9                   # Decimal hours between dumps of restart files.\n")
    f.write("    max_run_time:       70              # Maximal wall-clock time in hours. The application will exit when this limit is reached.\n")
    f.write("    resubmit_on_exit:   1               # Whether to run a command when exiting after the time limit has been reached.\n")
    f.write("    resubmit_command:   ./auto_resubmit.sh    # Command to run when time limit is reached. Compulsory if resubmit_on_exit is switched on. Potentially unsafe.\n\n")

    f.write("# Parameters governing domain decomposition\n")
    f.write("DomainDecomposition:\n")
    f.write("    trigger:        0.1                 # Fractional (<1) CPU time difference between MPI ranks required to trigger a new decomposition, or number of steps (>1) between decompositions\n\n")


    f.write("# Parameters for the hydrodynamics scheme\n")
    f.write("SPH:\n")
    f.write("    resolution_eta:   " + "{:.4f}".format(eta) +  "       # Target smoothing length in units of the mean inter-particle separation (1.2348 == 48Ngbs with the cubic spline kernel).\n")
    f.write("    delta_neighbours:   0.1             # The tolerance for the targetted number of neighbours.\n")
    f.write("    CFL_condition:      0.2             # Courant-Friedrich-Levy condition for time integration.\n")
    f.write("    h_max:              " + "{:.2f}".format(h_max) +  "             # Maximal allowed smoothing length (in internal units).\n")
    f.write("    viscosity_alpha:    1.5             # Override for the initial value of the artificial viscosity.\n\n")

    f.write("# Parameters for the self-gravity scheme\n")
    f.write("Gravity:\n")
    f.write("    eta:                0.025           # Constant dimensionless multiplier for time integration.\n")
    f.write("    MAC:                geometric        # Choice of mulitpole acceptance criterion: 'adaptive' OR 'geometric'.\n")
    f.write("    theta_cr:           0.7             # Opening angle for the purely geometric criterion.\n")
    if N == 1e5:
        f.write("    max_physical_baryon_softening: 0.02   # Physical softening length (in internal units).\n\n")
    elif N == 1e6:
        f.write("    max_physical_baryon_softening: 0.01   # Physical softening length (in internal units).\n\n")
    elif N == 1e7:
        f.write("    max_physical_baryon_softening: 0.005   # Physical softening length (in internal units).\n\n")
    else:
        raise ValueError("You messed up genius!")

    f.write("# Parameters for the task scheduling\n")
    f.write("Scheduler:\n")
    f.write("    max_top_level_cells:    48          # Maximal number of top-level cells in any dimension.\n")
    f.write("    cell_split_size:        400         # Maximal number of particles per cell (400 is the default value).\n")
    f.write("    task_per_cell:          2.0\n\n")

    f.write("# Parameters related to the equation of state\n")
    f.write("EoS:\n")
    if mat == "Til":
        f.write("    planetary_use_Til:  1 \n")
        f.write("    planetary_use_Til_iron:  1 \n")
        f.write("    planetary_use_Til_granite:  1 \n")
        f.write("    planetary_use_Til_basalt:  1 \n")
    elif mat == "ANEOS":
        f.write("    planetary_use_ANEOS_forsterite:  1 \n")
        f.write("    planetary_ANEOS_forsterite_table_file:   /cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/Planetary/EoSTables/ANEOS_forsterite_S19.txt \n")
        f.write("    planetary_use_ANEOS_Fe85Si15:  1 \n")
        f.write("    planetary_ANEOS_Fe85Si15_table_file:     /cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/Planetary/EoSTables/ANEOS_Fe85Si15_S20.txt \n")
        f.write("    planetary_use_ANEOS_iron:  1 \n")
        f.write("    planetary_ANEOS_iron_table_file:         /cosma/home/dp004/dc-ruiz1/data/swiftsim/examples/Planetary/EoSTables/ANEOS_iron_S20.txt \n")
    else:
        raise ValueError("You messed up genius!")

    f.close()
    
def make_script_slurm(N, m, **kwargs):
    
    time = kwargs.get("time", "72:00:00")
    hydro = kwargs.get("hydro")
    kernel = kwargs.get("kernel")
    
    # check time format
    error_message = f"time must be a string of format xx:xx:xx (hours:min:sec)"
    assert isinstance(time, str), error_message
    assert len(time) == 8, error_message
    assert time[2] == ":" and time[5] == ":", error_message
    for i in (0,1,3,4,6,7):
        assert time[i].isnumeric(), error_message
    
    folder_name = get_dir_settling(N, m, **kwargs)
    fn = os.path.join(folder_name, "script_slurm.sh")
    
    with open(fn, 'w') as f:
        f.write("#!/bin/bash -l\n\n")
        f.write("#   Number of processors\n")
        f.write("#SBATCH -N 1\n")
        f.write("#SBATCH --ntasks 28\n")
        f.write("#SBATCH --tasks-per-node=28\n")
        f.write("#SBATCH --exclusive\n")
        f.write("#   Job name\n")
        f.write("#SBATCH -J SWIFT\n")
        f.write("#   Output files\n")
        f.write("#SBATCH -o ./swift_%J.out\n")
        f.write("#SBATCH -e ./swift_%J.err\n")
        f.write("#   Queue\n")
        f.write("#SBATCH -p cosma7\n")
        f.write("#   Project\n")
        f.write("#SBATCH -A dp004\n")
        f.write("#   Email info\n")
        f.write("#SBATCH --mail-type=NONE\n")
        f.write("#SBATCH --mail-user=sergio.ruiz-bonilla@durham.ac.uk\n")
        f.write("#   Wall clock time limit (limit 72:00:00)\n")
        f.write("#SBATCH -t " + time + "\n\n")

        f.write("#   Ensure that the right modules are loaded\n")
        f.write("module purge\n")
        f.write("module load intel_comp/2020-update2 \n")
        f.write("module load intel_mpi/2020-update2\n")
        f.write("module load ucx/1.8.1 \n")
        f.write("module load parmetis/4.0.3-64bit \n")
        f.write("module load parallel_hdf5/1.10.6 \n")
        f.write("module load fftw/3.3.8cosma7    # On cosma 5 or 6 use fftw/3.3.8, on cosma 7 use fftw/3.3.8cosma7\n")
        f.write("module load gsl/2.5 \n")
        f.write("module load llvm/10.0.1         # Only necessary if wanting to use the code formatting tool\n\n")

        swift = get_swift_flavour(hydro, kernel)
        f.write(swift + " -G -s -a -t 28 planet.yml\n\n")
        f.write("#   Print info\n")
        f.write("echo \"Job done! \"\n")
        f.write("sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode\n")
        f.write("exit\n")
        
    fn = os.path.join(folder_name, "script_slurm_restart.sh")
    
    with open(fn, 'w') as f:
        f.write("#!/bin/bash -l\n\n")
        f.write("#   Number of processors\n")
        f.write("#SBATCH -N 1\n")
        f.write("#SBATCH --ntasks 28\n")
        f.write("#SBATCH --tasks-per-node=28\n")
        f.write("#SBATCH --exclusive\n")
        f.write("#   Job name\n")
        f.write("#SBATCH -J SWIFT\n")
        f.write("#   Output files\n")
        f.write("#SBATCH -o ./swift_%J.out\n")
        f.write("#SBATCH -e ./swift_%J.err\n")
        f.write("#   Queue\n")
        f.write("#SBATCH -p cosma7\n")
        f.write("#   Project\n")
        f.write("#SBATCH -A dp004\n")
        f.write("#   Email info\n")
        f.write("#SBATCH --mail-type=NONE\n")
        f.write("#SBATCH --mail-user=sergio.ruiz-bonilla@durham.ac.uk\n")
        f.write("#   Wall clock time limit (limit 72:00:00)\n")
        f.write("#SBATCH -t 72:00:00\n\n")

        f.write("#   Ensure that the right modules are loaded\n")
        f.write("module purge\n")
        f.write("module load intel_comp/2020-update2 \n")
        f.write("module load intel_mpi/2020-update2\n")
        f.write("module load ucx/1.8.1 \n")
        f.write("module load parmetis/4.0.3-64bit \n")
        f.write("module load parallel_hdf5/1.10.6 \n")
        f.write("module load fftw/3.3.8cosma7    # On cosma 5 or 6 use fftw/3.3.8, on cosma 7 use fftw/3.3.8cosma7\n")
        f.write("module load gsl/2.5 \n")
        f.write("module load llvm/10.0.1         # Only necessary if wanting to use the code formatting tool\n\n")

        swift = get_swift_flavour(hydro, kernel)
        f.write(swift + " -G -s -a -t 28 -r planet.yml\n\n")
        f.write("#   Print info\n")
        f.write("echo \"Job done! \"\n")
        f.write("sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,MaxRSS,Elapsed,ExitCode\n")
        f.write("exit\n")
        
    return None

def start_settling_sim(N, m, **kwargs):
    # why doesn't work????
    
    # go to folder
    command = "cd " + get_dir_settling(N, m, **kwargs)
    print(command)
    #os.system(command)
    
    # give execution permision
    #command = "chmod u+x script_slurm.sh"
    #print(command)
    #os.system(command)
    
    # submit job
    command = "sbatch script_slurm.sh"
    print(command)
    #os.system(command)