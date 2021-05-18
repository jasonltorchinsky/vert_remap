#!/bin/bash

#SBATCH --time=0-0:00:30
#SBATCH --job-name=vert_remap
#SBATCH --output=%x.out
#SBATCH --error=%x.err

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
module load cray-netcdf

build_dir=../build

echo '-- Running vertical remap code...'
$build_dir/vert_remap ncell 4

echo '-- Vertical remap code complete!'
