#!/bin/bash
#SBATCH -J trail_1     # name of the job
#SBATCH -p shared # name of the partition: available options "standard, standard-low, gpu, hm"
#SBATCH --nodes=1
#SBATCH -n 16
#SBATCH -t 8:00:00    # walltime in HH:MM:SS, Max value 72:00:00
#SBATCH --mail-user=yourmail@mail.com
#SBATCH --mail-type=ALL


#list of modules you want to use, for example
#module load apps/python-package/python/conda-python/3.7_new
module load orca/421
module load xtb/650
source /home/${USER}/.bashrc
conda activate base

export PYTHONPATH=/home/18cy91r30/autore2/:$PYTHONPATH
export PYTHONPATH=/home/18cy91r30/autore2/interface:$PYTHONPATH
echo $PYTHONPATH >> env.log
which python3 >> env.log

python3 -u ~/rksn_automate/start.py -m start.xyz -r act_atom.txt > output.log
bash path.txt


