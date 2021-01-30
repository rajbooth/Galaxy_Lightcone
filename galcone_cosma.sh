#!/bin/bash -l

#SBATCH --ntasks 1
#SBATCH -J galcone
#SBATCH -o output/galcone.%J.out
#SBATCH -p cosma6
#SBATCH -A dp004
#SBATCH --mem-per-cpu=4096
#SBATCH -t 24:00:00
#SBATCH --mail-type=END                          # notifications for job done & fail
#SBATCH --mail-user=rb460@sussex.ac.uk

module purge
#load the modules used to build your program.
module load pythonconda3/4.5.4 

# Run the program
python make_galaxy_lightcone.py

