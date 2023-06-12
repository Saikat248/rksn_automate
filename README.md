# rksn_automate

This repository contains automation code developed for the research paper titled "Exploring Enantioselectivity in Organocatalysis: Insights from Bonding and Energy Decomposition Analyses". 


The code automates the TS search method for all the conformers based on the product structure.

# Dependencies

- Python (version >= 3.7 should work fine)
- No python external library required.
- Mainly tested with Anaconda Python.
- Without Anaconda the automation script should run fine but need to change couple of things (see below).

## External Programmes

- Interfaced with ORCA 4.2.1
- XTB and CREST from Grimme Lab 


# Installation

- Clone the `rksn_automate` repository.

- Set the absolute path in [lib.py](https://github.com/Saikat248/rksn_automate/blob/master/interface/lib.py) `def slurm_input(inp_file, jobname)` function.

- This code runs in the SLURM  queuing system. A slurm submit script `submit.sh` is provided. Please change the file accordingly. 

- If other than SLURM queuing system 
is used change the [lib.py](https://github.com/Saikat248/rksn_automate/blob/master/interface/lib.py) `def slurm_input(inp_file, jobname)` function and modify the `submit.sh` file. If you are not using Anaconda python, find the proper python interpreter path and modify this two file.

- In our test system both ORCA and XTB software are set up with the Module system.
Load proper version of ORCA and XTB in the `submit.sh` file. Without the module 
system also it is easy to set up the automation script. 

# Running the Script 

The run command is provided inside the `submit.sh` file.

```bash
python3 -u ~/rksn_automate/start.py -m start.xyz -r act_atom.txt > output.log
```
where `rksn_automate` repo is cloned in the $HOME directory of the user.


The run command required a `act_atom.txt` file and an optional `constrain` file.
Samples for the files are provided in the repository.

The `act_atom.txt` file contains the index of the two atoms where bond breaking or forming will happen according to the `start.xyz` file.
The `constrain` file contains the the Dihedral Constrain if required during the bond scan 
between the two active atoms. This format of the file is according to the ORCA input file.
All the index starts from 0.

The `start.py` file first do the conformation search using CREST and XTB. This script will 
generate a `path.txt` file which will be submitted later automatically for parallel run of the 
all conformers. 

The DFT methods and settings for ORCA runs are provided at the beginning of the [lib.py](https://github.com/Saikat248/rksn_automate/blob/master/interface/lib.py) file. The SMD single point methods are also present. You can set up your preferred DFT Functional and basis set here.
Change `PROCS` variable for setting up the number of ORCA parallel processors. The `BLOCKS` variable determines how many low energy conformers will be taken from the crest conformer search file for the reaction path study.


# Contact

For any inquiries or questions regarding the code or research paper, please contact saikat403@gamil.com.






