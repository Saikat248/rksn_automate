import os
import glob
import shutil
from interface.lib import * 
import time

"""
Need to export this PYTHONPATH before running 
export PYTHONPATH="full-directory-path/rksn_automate":$PYTHONPATH
"""

cpu_st = time.process_time()
wall_st = time.time()

cwd = os.getcwd()

# Finding Active Atom from rule.txt file 

with open('act_atom.txt', 'r') as fp:
    lines = fp.readlines()

active_atom = []
for line in lines:
    for i in (line.split(',')):
        active_atom.append(int(i))

pre_complex = glob.glob('conf*.xyz')

# Check for constrain

isConstrain =  os.path.isfile('constrain')


# For Multiple ORCA constrain handeling
if isConstrain:
    with open('constrain') as fp:
        lines = fp.readlines()
    constrain = lines
else:
    print('Constrain File Not Found')
    constrain = None

# Doing Initial Crude Optimization

os.makedirs(cwd+'/initial_precomplex_opt')
initial_opt_path = cwd+'/initial_precomplex_opt'
shutil.copy(cwd+'/'+pre_complex[0], initial_opt_path)
os.chdir(initial_opt_path)
orca_opt_inp(pre_complex[0])
orca_run('opt.inp', 'result.out', 'opt_'+pre_complex[0])
os.chdir(cwd)

# Doing SCAN 
os.makedirs(cwd+'/scan')
scan_path = cwd+'/scan'
shutil.copy(initial_opt_path+'/opt.xyz', scan_path)
os.chdir(scan_path)
# orca_scan_inp('opt.xyz', active_atom[0], active_atom[1])
orca_scan_constrain_inp('opt.xyz', active_atom[0], active_atom[1], constrain=constrain)
orca_run('scan.inp', 'result.out', 'scan_'+pre_complex[0])
os.chdir(cwd)

# Doing NEB
os.makedirs(cwd+'/neb')
neb_path = (cwd+'/neb')
shutil.copy(initial_opt_path+'/opt.xyz', neb_path)
shutil.copy(scan_path+'/scan.xyz', neb_path)
os.chdir(neb_path)
os.rename('opt.xyz', 'end.xyz')
os.rename('scan.xyz', 'start.xyz')
orca_neb_inp('start.xyz', 'end.xyz')
orca_run('neb.inp', 'result.out', 'neb_'+pre_complex[0])
os.chdir(cwd)

# Doing TS Opt and IRC

os.makedirs(cwd+'/ts_opt')
ts_path = cwd+'/ts_opt'
shutil.copy(neb_path+'/neb_TSOpt.xyz', ts_path)
os.chdir(ts_path)
orca_ts_irc_inp('{ '+str(active_atom[0])+' '+str(active_atom[1])+' }', 'neb_TSOpt.xyz')
orca_ts_irc_inp(f'{{{active_atom[0]} {active_atom[1]}}}', 'neb_TSOpt.xyz')
orca_run('ts.inp', 'result.out', 'TS Opt_IRC'+pre_complex[0])
os.makedirs(ts_path+'/start')
os.makedirs(ts_path+'/pdt')
for j in os.listdir():
    if j == 'ts_IRC_F.xyz' or j == 'ts_IRC_B.xyz':
        print(f'********** IRC End Points are found *************')
        dist = get_dist(j, active_atom[0], active_atom[1])
        if dist < 1.85:
            shutil.copy(ts_path+'/'+j, ts_path+'/pdt')
            os.chdir(ts_path+'/pdt')
            orca_opt_freq_inp(j)
            orca_run('opt.inp', 'result.out', 'pdt_'+pre_complex[0])
            os.chdir(ts_path)
        else:
            shutil.copy(ts_path+'/'+j, ts_path+'/start')
            os.chdir(ts_path+'/start')
            orca_opt_freq_inp(j)
            orca_run('opt.inp', 'result.out', 'start_'+pre_complex[0])
            os.chdir(ts_path)

# Running SMD Solvation

"""
Comment the Whole Code (TS-SinglePoint -- Pdt-SinglePoint) if not
Required
"""

# TS-SinglePoint

os.chdir(ts_path)
os.makedirs(ts_path+'/ts_singlepoint')
ts_sing_path = ts_path+'/ts_singlepoint' 
shutil.copy('ts.xyz', ts_sing_path)
os.chdir(ts_sing_path)
orca_smd_inp('ts.xyz')
orca_run('smd.inp', 'result.out', 'ts_'+pre_complex[0])

# Start Singlepoint
os.chdir(ts_path+'/start')
os.makedirs(ts_path+'/start/start_singlepoint')
shutil.copy('opt.xyz', ts_path+'/start/start_singlepoint')
os.chdir(ts_path+'/start/start_singlepoint')
orca_smd_inp('opt.xyz')
orca_run('smd.inp', 'result.out', 'start_'+pre_complex[0])

# Pdt SinglePoint
os.chdir(ts_path+'/pdt')
os.makedirs(ts_path+'/pdt/pdt_singlepoint')
shutil.copy('opt.xyz', ts_path+'/pdt/pdt_singlepoint')
os.chdir(ts_path+'/pdt/pdt_singlepoint')
orca_smd_inp('opt.xyz')
orca_run('smd.inp', 'result.out', 'pdt_'+pre_complex[0])

cpu_et = time.process_time()
wall_et = time.time()

wall_elapsed = wall_et - wall_st
cpu_elapsed = cpu_et - cpu_st

print(2*'\n')
print('Total Timings: '+'\n')
print('==========================')
print('Wall Time: ', time.strftime("%H:%M:%S", time.gmtime(wall_elapsed)))
print('CPU Time: ', time.strftime("%H:%M:%S", time.gmtime(cpu_elapsed)))
print('==========================')


