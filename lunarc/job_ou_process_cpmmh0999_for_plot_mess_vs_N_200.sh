#!/bin/bash

# This script creates the run-scripts for the 4 runs with different number of threads

# set number of threads
seed_data_1=100
seed_data_2=200
seed_data_3=300
seed_data_4=400
seed_data_5=500
# this will be extended later on...


# set FILE names
FILE1="job_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seed_data_1}.sh"
FILE2="job_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seed_data_2}.sh"
FILE3="job_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seed_data_3}.sh"
FILE4="job_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seed_data_4}.sh"
FILE5="job_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seed_data_5}.sh"


# arrays with file names and nbr of threads
FILES=($FILE1 $FILE2 $FILE3 $FILE4 $FILE5)
seeds_for_data=($seed_data_1 $seed_data_2 $seed_data_3 $seed_data_4 $seed_data_5)

# loop over each file and thread
for ((i=0;i<=$((${#FILES[@]} - 1));i++)); do

# check if file exists, if not create an empty file
#if [ ! -e ${FILES[$i]} ]; then
#  echo >> ${FILES[$i]}
#fi


# create empty file
echo >> ${FILES[$i]}

outputfile="lunarc_output/outputs_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seeds_for_data[$i]}_%j.out"
errorfile="lunarc_output/errors_ou_process_cpmmh0999_for_plot_mess_vs_N_200_${seeds_for_data[$i]}_%j.err"

cat > ${FILES[$i]} << EOF
#!/bin/sh


# Set up for run:

# need this since I use a LU project
# #SBATCH -A snic2019-3-630 #TODO this should be updated!

#SBATCH -A lu2019-2-19
#SBATCH -p lu

# use gpu nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive

# time consumption HH:MM:SS
#SBATCH -t 10:00:00

# name for script
#SBATCH -J ou_cpmmh_0999_200

# controll job outputs
#SBATCH -o lunarc_output/outputs_ou_cpmmh_%j.out
#SBATCH -e lunarc_output/errors_ou_cpmmh_%j.err

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
julia /home/samwiq/'SDEMEM and CPMMH'/SDEMEM_and_CPMMH/src/'SDEMEM OU process'/run_script_cpmmh_for_plot_mess_vs_N.jl 200 0.999 ${seeds_for_data[$i]} # M_mixtures N_time nbr_particles correlation seed
EOF


# run job
sbatch ${FILES[$i]}


done
