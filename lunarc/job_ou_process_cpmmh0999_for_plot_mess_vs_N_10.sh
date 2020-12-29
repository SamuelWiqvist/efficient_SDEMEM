#!/bin/bash


# set number of threads
seed_data_1=100
seed_data_2=200
seed_data_3=300
seed_data_4=400
seed_data_5=500
seed_data_6=600
seed_data_7=700
seed_data_8=800
seed_data_9=900
seed_data_10=10
seed_data_11=11
seed_data_12=12
seed_data_13=13
seed_data_14=14
seed_data_15=15
seed_data_16=16
seed_data_17=17
seed_data_18=18
seed_data_19=19
seed_data_20=20
seed_data_21=21
seed_data_22=22
seed_data_23=23
seed_data_24=24
seed_data_25=25

# set FILE names
FILE1="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_1}.sh"
FILE2="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_2}.sh"
FILE3="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_3}.sh"
FILE4="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_4}.sh"
FILE5="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_5}.sh"
FILE6="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_6}.sh"
FILE7="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_7}.sh"
FILE8="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_8}.sh"
FILE9="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_9}.sh"
FILE10="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_10}.sh"
FILE11="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_11}.sh"
FILE12="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_12}.sh"
FILE13="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_13}.sh"
FILE14="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_14}.sh"
FILE15="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_15}.sh"
FILE16="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_16}.sh"
FILE17="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_17}.sh"
FILE18="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_18}.sh"
FILE19="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_19}.sh"
FILE20="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_20}.sh"
FILE21="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_21}.sh"
FILE22="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_22}.sh"
FILE23="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_23}.sh"
FILE24="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_24}.sh"
FILE25="job_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seed_data_25}.sh"

# arrays with file names and nbr of threads
FILES=($FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7 $FILE8 $FILE9 $FILE10 $FILE11 $FILE12 $FILE13 $FILE14 $FILE15 $FILE16 $FILE17 $FILE18 $FILE19 $FILE10 $FILE21 $FILE22 $FILE23 $FILE24 $FILE25) 
seeds_for_data=($seed_data_1 $seed_data_2 $seed_data_3 $seed_data_4 $seed_data_5 $seed_data_6 $seed_data_7 $seed_data_8 $seed_data_9 $seed_data_10 $seed_data_11 $seed_data_12 $seed_data_13 $seed_data_14 $seed_data_15 $seed_data_16 $seed_data_17 $seed_data_18 $seed_data_19 $seed_data_20 $seed_data_21 $seed_data_22 $seed_data_23 $seed_data_24 $seed_data_25) 

# loop over each file and thread
for ((i=0;i<=$((${#FILES[@]} - 1));i++)); do

# check if file exists, if not create an empty file
#if [ ! -e ${FILES[$i]} ]; then
#  echo >> ${FILES[$i]}
#fi


# create empty file
echo >> ${FILES[$i]}

outputfile="lunarc_output/outputs_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seeds_for_data[$i]}_%j.out"
errorfile="lunarc_output/errors_ou_process_cpmmh0999_for_plot_mess_vs_N_10_${seeds_for_data[$i]}_%j.err"

cat > ${FILES[$i]} << EOF
#!/bin/sh


# Set up for run:

# need this since I use a LU project
# #SBATCH -A lu2020-2-7
# #SBATCH -p lu

#SBATCH -A snic2019-3-630 ##SBATCH -A snic2016-x-xxx


# use gpu nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive

# time consumption HH:MM:SS
#SBATCH -t 10:00:00

# name for script
#SBATCH -J ou_cpmmh_0999_10

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
julia /home/samwiq/'SDEMEM and CPMMH'/SDEMEM_and_CPMMH/src/'SDEMEM OU process'/run_script_cpmmh_for_plot_mess_vs_N.jl 10 0.999 ${seeds_for_data[$i]} # M_mixtures N_time nbr_particles correlation seed
EOF


# run job
sbatch ${FILES[$i]}


done
