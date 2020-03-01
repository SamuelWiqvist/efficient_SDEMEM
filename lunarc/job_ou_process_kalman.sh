#!/bin/sh


# Set up for run:

# need this since I use a LU project
#SBATCH -A snic2019-3-630


# use gpu nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive

# time consumption HH:MM:SS
#SBATCH -t 01:00:00

# name for script
#SBATCH -J ou_kalman

# controll job outputs
#SBATCH -o lunarc_output/outputs_ou_kalman_%j.out
#SBATCH -e lunarc_output/errors_ou_kalman_%j.err

# notification
#SBATCH --mail-user=samuel.wiqvist@matstat.lu.se
#SBATCH --mail-type=ALL

# load modules

ml load GCC/6.4.0-2.28
ml load OpenMPI/2.1.2
ml load julia/1.0.0

# set correct path
pwd
cd ..
pwd

export JULIA_NUM_THREADS=1

# run program
julia /home/samwiq/'SDEMEM and CPMMH'/SDEMEM_and_CPMMH/src/'SDEMEM OU process'/run_script_kalman.jl $1
