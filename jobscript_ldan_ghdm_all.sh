#!/bin/bash

#This template is adapted from Molly's foggie/examples/RunScript.sh; see there for explanations

#PBS -N ldan_ghdm_all

##### Resources #####
#PBS -W group_list=s2358
#PBS -l select=1:ncpus=16:mpiprocs=16:mem=1500GB
#PBS -l walltime=72:00:00

##### Queue #####
#PBS -q ldan

##### Mail Options #####
#PBS -m abe
#PBS -M ayan.acharyya@inaf.it

#set output and error directories
#PBS -j oe
#PBS -o pbs_ldan_ghdm_all.out

###### Load necessary modules #########
module purge
module load mpi-hpe/mpt.2.30
source /home5/aachary2/miniconda3/etc/profile.d/conda.sh
conda activate py380

export PYTHONPATH=$PYTHONPATH:/home5/aachary2/miniconda3/envs/py380/bin/

# ------------------------

###### Load necessary modules #########
module purge
module load mpi-hpe/mpt.2.30
source /home5/aachary2/miniconda3/etc/profile.d/conda.sh
conda activate py380

export PYTHONPATH=$PYTHONPATH:/home5/aachary2/miniconda3/envs/py380/bin/

# ------------------------

##### Change to current working directory #####
cd /nobackupp19/aachary2/foggie_craft/pleiades_workdir

##### Execute Program #####
top -b -u $USER -n 259200 1>top_ldan_ghdm_all.out &
mpiexec -n 16 python /nobackupp19/aachary2/foggie_craft/foggie_craft/gethostdm_all.py losdm 100 200 0.5 1>output_ldan_ghdm_all.out 2>&1
killall 'top'