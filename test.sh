#!/bin/bash -l

#SBATCH --ntasks 1
#SBATCH -J lightcone
#SBATCH -o output/lightcone.%J.out
#SBATCH -p cosma6
#SBATCH -A dp004
#SBATCH --mem-per-cpu=64800
#SBATCH -t 24:00:00
#SBATCH --mail-type=END                          # notifications for job done & fail
#SBATCH --mail-user=rb460@sussex.ac.uk

module purge

source activate python3_env

#load the modules used to build your program.
module load intel_comp/2018
module load intel_mpi/2018

# Run the program
python batch_test.py

