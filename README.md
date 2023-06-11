# rksn_automate

This repository contains automation code developed for the research paper titled "Exploring Enantioselectivity in Organocatalysis: Insights from Bonding and Energy Decomposition Analyses". 


The code automates the TS search method for all the conformers based on the product structure.

# Dependencies

Python (version >= 3.7 should work fine)
No python external library required.
Mainly tested with Anaconda Python.
Without Anaconda the automation script should run fine but need to change couple of things (see below).

## External Programmes

- Interfaced with ORCA 4.2.1
- XTB and CREST from Grimme Lab 


# Installation

Clone the `rksn_automate` repository.
Then set the absolute path in [https://github.com/Saikat248/rksn_automate/blob/master/interface/lib.py](lib.py) `def slurm_input(inp_file, jobname)` function.
This code runs in the SLURM  queuing system.
A slurm submit script `submit.sh` is provided. Please change the file accordingly. 


# Contact

For any inquiries or questions regarding the code or research paper, please contact saikat403@gamil.com.






