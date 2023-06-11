import glob
import os
import shutil
import argparse
import time
from interface.lib import *


# cwd = os.getcwd()
# pdt_complex = 'start.xyz'

# XTB Conformer search

def initial_run(pdt_complex, rule_file, constrain_file=None):
    cwd = os.getcwd()
    crest_run(pdt_complex, 'result.out')
    crest_conf_split()
    os.makedirs(cwd+'/conformations')
    conf_lst = glob.glob('conf*.xyz')
    print(conf_lst)
    isConstrain =  os.path.isfile('constrain')
    paths = []
    for i in conf_lst:
        os.makedirs(cwd+'/conformations/'+(i.split('.')[0]))
        shutil.copy(i, cwd+'/conformations/'+(i.split('.')[0]))
        shutil.copy('/home/18cy91r30/autore2/driver.py', cwd+'/conformations/'+(i.split('.')[0]))
        shutil.copy(cwd+'/'+rule_file, cwd+'/conformations/'+(i.split('.')[0]))
        if isConstrain:
            shutil.copy(cwd+'/constrain', cwd+'/conformations/'+(i.split('.')[0]))
        else:
            pass
        os.chdir(cwd+'/conformations/'+(i.split('.')[0]))
        slurm_input('submit.sh',i.split('.')[0]) 
        path = os.getcwd()
        paths.append(path) 
        os.chdir(cwd)

    # print(paths)

    with open('path.txt', 'w') as fp:
        for i in paths:
            fp.writelines('cd '+i+'\n')
            fp.writelines('sbatch submit.sh'+'\n')




def main():
    parser = argparse.ArgumentParser(description='Parser for autore2')
    parser.add_argument('-m', '--molecule', type=str, nargs=1, required=True, help='Name of the pdt')
    parser.add_argument('-r', '--rule', type=str, nargs=1, required=True, help='Active atom file name has to be act_atom.txt')

    args = parser.parse_args()
    pdt_complex = args.molecule[0]
    rule_file = args.rule[0]

    cpu_st = time.process_time()
    wall_st = time.time()

    initial_run(pdt_complex, rule_file)
    # print(pdt_complex)
    # print(rule_file)

    cpu_et = time.process_time()
    wall_et = time.time()

    wall_elapsed = wall_et - wall_st
    cpu_elapsed = cpu_et - cpu_st
    
    print(2*'\n')
    print('Total Timings: '+'\n')
    print('==========================')
    print('Wall Time: ', time.strftime("%H:%M:%S", time.gmtime(wall_elapsed)))
    print('CPU Time:  ', time.strftime("%H:%M:%S", time.gmtime(cpu_elapsed)))
    print('==========================')

if __name__ == '__main__':
    main()
