#!/bin/bash

#SBATCH -JNAME
#SBATCH --account=hive-gli323
#SBATCH -t30:00:00
#SBATCH -phive
#SBATCH -ooutput.out

cd WORKDIR
module load anaconda3/2022.05
conda activate /storage/home/hhive1/jbrandt35/.conda/envs/rebound_gr
srun python Simulation.py