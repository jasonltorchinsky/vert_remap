#!/bin/bash

#SBATCH --time=0-0:01:00
#SBATCH --job-name=vert_remap
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --tasks-per-node=01
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jason.torchinsky@wisc.edu
#SBATCH --mail-type=none

# Set up for run

echo '-- Setting submission variables...'
cd $SLURM_SUBMIT_DIR
module purge
module load PrgEnv-intel

build_dir=../build

echo '-- Running vertical remap code...'
$build_dir/vert_remap nlev 1001

