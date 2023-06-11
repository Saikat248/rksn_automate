import os
import sys
import subprocess as subp
import numpy as np


# Default Parameters for ORCA and CREST Calculations

DFT_METHOD = 'B97-D3'
BASIS_SET = 'def2-SVP def2-SVP/C def2/J'
H_BASIS_SET = 'def2-SVP def2-SVP/C def2/J'
CHARGE = '0'
SPIN = '1'
PROCS = '16'
O_PARAMS = 'RIJCOSX Grid6 NormalSCF NoPop NoFinalGrid'
BLOCKS = 24
DRIVER_TIME='70:00:00'


# Single Point ORCA Setup

SMD_DFT_METHOD = 'M06'
SMD_BASIS_SET = 'def2-TZVP def2-TZVP/C def2/J'
SMD_O_PARAMS = 'RIJCOSX Grid6 NormalSCF NoPop NoFinalGrid'
SMD_SOLVENT = 'METHANOL'

# Get Full Path of an system executable 

def which(program):
    """
    Find the path of an executable program.

    Parameters:
        program (str): The name of the program or executable file.

    Returns:
        str: The full path to the executable if found, None otherwise.
    """
    def is_exe(exec_path):
        """
        Check if a given file path is an executable.

        Parameters:
            exec_path (str): The file path to check.

        Returns:
            bool: True if the file is executable, False otherwise.
        """
        return os.path.isfile(exec_path) and os.access(exec_path, os.X_OK)

    file_path, file_name = os.path.split(program)
    if file_path:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

def orca_neb_inp(start_xyz, end_xyz):
    """
    Generate an ORCA input file for performing a NEB (Nudged Elastic Band) calculation.

    Parameters:
        start_xyz (str): File path to the starting structure in XYZ format.
        end_xyz (str): File path to the ending structure in XYZ format.

    Writes:
        Creates a new file 'neb.inp' with the generated ORCA input.

    Returns:
        None
    """
    with open('neb.inp', 'w') as nebfile:
        nebfile.writelines(f'!NEB-TS Opt {DFT_METHOD} {BASIS_SET} {O_PARAMS}' + '\n')
        nebfile.writelines('%scf' + '\n')
        nebfile.writelines('   MaxIter 999' + '\n')
        nebfile.writelines('   CNVDIIS true' + '\n')
        nebfile.writelines('   CNVSOSCF true' + '\n')
        nebfile.writelines('   AutoStart true' + '\n')
        nebfile.writelines('end' + '\n')
        nebfile.writelines('%neb' + '\n')
        nebfile.writelines(f'   neb_end_xyzfile "{end_xyz}"'+'\n')
        nebfile.writelines('   NImages 14'+'\n')
        nebfile.writelines('   PrintLevel 1'+'\n')
        nebfile.writelines('   PreOpt_Ends False'+'\n')
        nebfile.writelines('   MaxIter 1000'+'\n')
        nebfile.writelines('end'+'\n')
        nebfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        nebfile.writelines(f'* xyzfile {CHARGE} {SPIN} {start_xyz}'+'\n')

def orca_opt_freq_inp(xyz_file):
    """
    Generate an ORCA input file for performing geometry optimization and frequency calculation.

    Parameters:
        xyz_file (str): File path to the structure in XYZ format.

    Writes:
        Creates a new file 'opt.inp' with the generated ORCA input.

    Returns:
        None
    """
    with open('opt.inp', 'w') as optfile:
        optfile.writelines(f'!Opt NumFreq {DFT_METHOD} {H_BASIS_SET} {O_PARAMS}' + '\n')
        optfile.writelines('%scf' + '\n')
        optfile.writelines('   MaxIter 999' + '\n')
        optfile.writelines('   CNVDIIS true' + '\n')
        optfile.writelines('   CNVSOSCF true' + '\n')
        optfile.writelines('   AutoStart true' + '\n')
        optfile.writelines('end' + '\n')
        optfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        optfile.writelines(f'* xyzfile {CHARGE} {SPIN} {xyz_file}'+'\n')

def orca_opt_inp(xyz_file):
    """
    Generate an ORCA input file for performing geometry optimization.

    Parameters:
        xyz_file (str): File path to the structure in XYZ format.

    Writes:
        Creates a new file 'opt.inp' with the generated ORCA input.

    Returns:
        None
    """

    with open('opt.inp', 'w') as optfile:
        optfile.writelines(f'!Opt {DFT_METHOD} {BASIS_SET} {O_PARAMS}' + '\n')
        optfile.writelines('%scf' + '\n')
        optfile.writelines('   MaxIter 999' + '\n')
        optfile.writelines('   CNVDIIS true' + '\n')
        optfile.writelines('   CNVSOSCF true' + '\n')
        optfile.writelines('   AutoStart true' + '\n')
        optfile.writelines('end' + '\n')
        optfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        optfile.writelines(f'* xyzfile {CHARGE} {SPIN} {xyz_file}'+'\n')



def xtbscan_split(scan_log_file):
    """
    Split an XTB Scan log file into separate start and end XYZ files.

    Parameters:
        scan_log_file (str): File path to the XTBScan log file.

    Writes:
        Creates two new files: 'end.xyz' and 'start.xyz' with the extracted data.

    Returns:
        None
    """
    with open(scan_log_file, 'r') as fp:
        l = fp.readlines()
        number_of_atoms = int(l[0])
    with open('end.xyz', 'w') as endfile:
        for i in l[0:number_of_atoms+2]:
            endfile.writelines(i)
    with open('start.xyz', 'w') as startfile:
        for i in (reversed(l[::-1][0:number_of_atoms+2])):
            startfile.writelines(i)


def orca_ts_irc_inp(active_atoms_dict, xyz_file):
    """
    Generate an ORCA input file for performing Transition State (TS) optimization and IRC (Intrinsic Reaction Coordinate) calculation.

    Parameters:
        active_atoms_dict (str): A dictionary specifying the active atoms for TS optimization and IRC.
        xyz_file (str): File path to the structure in XYZ format.

    Writes:
        Creates a new file 'ts.inp' with the generated ORCA input.

    Returns:
        None
    """
    with open('ts.inp', 'w') as tsfile:
        tsfile.writelines(f'! OptTS NumFreq {DFT_METHOD} {H_BASIS_SET} {O_PARAMS}' + '\n')
        tsfile.writelines('%scf' + '\n')
        tsfile.writelines(' MaxIter 999' + '\n')
        tsfile.writelines('   CNVDIIS true' + '\n')
        tsfile.writelines('   CNVSOSCF true' + '\n')
        tsfile.writelines('   AutoStart true' + '\n')
        tsfile.writelines('end' + '\n')
        tsfile.writelines('%geom' + '\n')
        tsfile.writelines('   TS_Active_Atoms ' + active_atoms_dict + '\n')
        tsfile.writelines('     end'+'\n')
        tsfile.writelines('   Calc_Hess true'+'\n')
        tsfile.writelines('   Hybrid_Hess ' +
                          active_atoms_dict + ' end' + '\n')
        tsfile.writelines('end'+'\n' + '\n')
        tsfile.writelines('%irc'+'\n')
        tsfile.writelines('   MaxIter 15'+'\n')
        tsfile.writelines('   PrintLevel 1'+'\n')
        tsfile.writelines('   Direction both'+'\n')
        tsfile.writelines('end'+'\n' + '\n')
        tsfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        tsfile.writelines(f'* xyzfile {CHARGE} {SPIN} {xyz_file}'+'\n')


def get_dist(xyz_file, atom1, atom2):
    """
    Calculate the Euclidean distance between two atoms in an XYZ file.

    Parameters:
        xyz_file (str): File path to the structure in XYZ format.
        atom1 (int): Index of the first atom.
        atom2 (int): Index of the second atom.

    Returns:
        float: The Euclidean distance between the two atoms.
    """
    with open(xyz_file, 'r') as xyzfile:
        lines = xyzfile.readlines()[2:]
    fp = lines[atom1].split()[1:]
    sp = lines[atom2].split()[1:]
    fp = np.array([float(i) for i in fp])
    sp = np.array([float(j) for j in sp])
    return np.linalg.norm(fp - sp)


def crest_run(xyzfile, output):
    """
    Run CREST with the specified XYZ file as input and save the output to a file.

    Parameters:
        xyzfile (str): File path to the input structure in XYZ format.
        output (str): File path to save the CREST output.

    Returns:
        None
    """
    crest_path = which('crest')
    # crest xyzfile -opt vtight -gfn2 -T procs
    with open(output, 'w') as ouputfile:
        run_crest = subp.Popen([crest_path, xyzfile, '-opt', 'vtight', '-gfn2', '-T', PROCS], stdout=ouputfile)
    run_crest.communicate()
    run_crest.poll()
    exitcode = run_crest.returncode
    if exitcode == 0:
        print(f'Crest Run Finished Successfully with exit code {exitcode}')
    else:
        print(f'*** WARNING! *** Crest Run Failure with exit code {exitcode}')
        sys.exit()


def crest_conf_split():
    """
    Split the CREST conformers XYZ file into separate XYZ files for each conformer.

    Parameters:
        None

    Writes:
        Creates multiple new files: 'conf1.xyz', 'conf2.xyz', ..., 'confN.xyz' with the splitted conformers.

    Returns:
        None
    """
    with open('crest_conformers.xyz', 'r') as fp:
        l = fp.readlines()
        number_of_atoms = int(l[0])
        j = 1
        for i in range(0, len(l), number_of_atoms+2):
            with open('conf' + str(j) + '.xyz', 'w') as new_file:
                new_file.writelines(l[i:i+number_of_atoms+2])
            j = j + 1
            if j == BLOCKS:
                break
    print('Conformers are splitted')

def crest_conf_all_split():
    with open('crest_conformers.xyz', 'r') as fp:
        l = fp.readlines()
        number_of_atoms = int(l[0])
        j = 1
        for i in range(0, len(l), number_of_atoms+2):
            with open('conf' + str(j) + '.xyz', 'w') as new_file:
                new_file.writelines(l[i:i+number_of_atoms+2])
            j = j + 1
    print('Conformers are splitted')

def orca_scan_inp(xyz_file, atom1, atom2):
    """
    Generate an ORCA input file for performing a bond scan.

    Parameters:
        xyz_file (str): File path to the input structure in XYZ format.
        atom1 (int): Index of the first atom involved in the bond scan.
        atom2 (int): Index of the second atom involved in the bond scan.

    Returns:
        None
    """
    initial_distance = get_dist(xyz_file, atom1, atom2)
    with open('scan.inp', 'w') as scanfile:
        scanfile.writelines(f'!Opt {DFT_METHOD} {BASIS_SET} {O_PARAMS}' + '\n')
        scanfile.writelines('%scf' + '\n')
        scanfile.writelines(' MaxIter 999' + '\n')
        scanfile.writelines('   CNVDIIS true' + '\n')
        scanfile.writelines('   CNVSOSCF true' + '\n')
        scanfile.writelines('   AutoStart true' + '\n')
        scanfile.writelines('end' + '\n')
        scanfile.writelines(f'%geom Scan' + '\n')
        scanfile.writelines(f'        B {atom1} {atom2} = {initial_distance:.3f}, 3.30, 10'+'\n')
        scanfile.writelines('        end'+'\n')
        scanfile.writelines('      end'+'\n')
        scanfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        scanfile.writelines(f'* xyzfile {CHARGE} {SPIN} {xyz_file}'+'\n')

def orca_scan_constrain_inp(xyz_file, atom1, atom2, constrain=None):
    """
    Generate an ORCA input file for performing a constrained bond scan.

    Parameters:
        xyz_file (str): File path to the input structure in XYZ format.
        atom1 (int): Index of the first atom involved in the bond scan.
        atom2 (int): Index of the second atom involved in the bond scan.
        constrain (list, optional): List of constraint strings for constraining specific coordinates during the scan.
                                   Default is None.

    Returns:
        None
    """
    initial_distance = get_dist(xyz_file, atom1, atom2)
    with open('scan.inp', 'w') as scanfile:
        scanfile.writelines(f'!Opt {DFT_METHOD} {BASIS_SET} {O_PARAMS}' + '\n')
        scanfile.writelines('%scf' + '\n')
        scanfile.writelines(' MaxIter 999' + '\n')
        scanfile.writelines('   CNVDIIS true' + '\n')
        scanfile.writelines('   CNVSOSCF true' + '\n')
        scanfile.writelines('   AutoStart true' + '\n')
        scanfile.writelines('end' + '\n')
        scanfile.writelines(f'%geom Scan' + '\n')
        scanfile.writelines(f'        B {atom1} {atom2} = {initial_distance:.3f}, 3.98, 10'+'\n')
        scanfile.writelines('        end'+'\n')
        if constrain is not None:
            scanfile.writelines(f'  Constraints' + '\n')
            # scanfile.writelines(f'  '+ constrain+'\n')
            for i in constrain:
                scanfile.writelines(f'  '+ i)
            scanfile.writelines('\n')
            scanfile.writelines(f'  end' + '\n')
        scanfile.writelines('end'+'\n')
        scanfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        scanfile.writelines(f'* xyzfile {CHARGE} {SPIN} {xyz_file}'+'\n')

def orca_smd_inp(xyz_file):
    """
    Generate an ORCA input file for performing a calculation with SMD solvation model.

    Parameters:
        xyz_file (str): File path to the input structure in XYZ format.

    Returns:
        None
    """
    with open('smd.inp', 'w') as optfile:
        optfile.writelines(f'!{SMD_DFT_METHOD} {SMD_BASIS_SET} {SMD_O_PARAMS}' + '\n')
        optfile.writelines('!CPCM' + '\n')
        optfile.writelines('%cpcm' + '\n')
        optfile.writelines('    smd true' + '\n')
        optfile.writelines(f'    SMDsolvent "{SMD_SOLVENT}"' + '\n')
        optfile.writelines('end' + '\n')
        optfile.writelines('%scf' + '\n')
        optfile.writelines('   MaxIter 999' + '\n')
        optfile.writelines('   CNVDIIS true' + '\n')
        optfile.writelines('   CNVSOSCF true' + '\n')
        optfile.writelines('   AutoStart true' + '\n')
        optfile.writelines('end' + '\n')
        optfile.writelines(f'%pal nprocs {PROCS} end'+'\n' + '\n')
        optfile.writelines(f'* xyzfile {CHARGE} {SPIN} {xyz_file}'+'\n')



def orca_run(inp_file, output, jobname):
    """
    Run an ORCA job using the specified input file.

    Parameters:
        inp_file (str): File path to the ORCA input file.
        output (str): File path to the output file for storing the ORCA job output.
        jobname (str): Name of the ORCA job.

    Returns:
        None
    """
    with open(output, 'w') as orca_output:
        run_job = subp.Popen([which('orca'), inp_file], stdout=orca_output)
    run_job.communicate()
    run_job.poll()
    exit_code = run_job.returncode
    if exit_code == 0:
        print(f'************ {jobname} Finished **************')
    else:
        print(f'************ {jobname} Failed with {exit_code} **************')
        sys.exit()



def slurm_input(inp_file, jobname):
    """
    Generate a SLURM input script for submitting a job.

    Parameters:
        inp_file (str): File path to the input file required by the job.
        jobname (str): Name of the job.

    Returns:
        None
    """
    with open(inp_file, 'w') as slurm_file:
        slurm_file.writelines('#!/bin/bash' + '\n')
        slurm_file.writelines(f'#SBATCH -J {jobname}'+'\n')
        slurm_file.writelines(f'#SBATCH -p shared'+'\n')
        slurm_file.writelines(f'#SBATCH --nodes=1'+'\n')
        slurm_file.writelines(f'#SBATCH -n {PROCS}'+'\n')
        slurm_file.writelines(f'#SBATCH -t {DRIVER_TIME}'+'\n')
        slurm_file.writelines(f'#SBATCH --mail-user=yourmail@mail.com'+'\n')
        slurm_file.writelines(f'#SBATCH --mail-type=ALL'+ 2*'\n')
        slurm_file.writelines(f'')
        slurm_file.writelines(f'module load orca/421'+ '\n')
        # Load Required python Conda Venv 
        slurm_file.writelines('source /home/${USER}/.bashrc'+'\n')
        slurm_file.writelines('export PYTHONPATH="/Full_path_to_rksn_automate":$PYTHONPATH'+'\n')
        slurm_file.writelines('export PYTHONPATH="/Full_path_to_rksn_automate/rksn_automate/interface":$PYTHONPATH'+'\n')
        slurm_file.writelines('conda activate base' +'\n')
        slurm_file.writelines(f'exe=$ORCAPATH/orca'+ '\n')
        slurm_file.writelines(f'python3 -u driver.py > driver.log'+ '\n')

