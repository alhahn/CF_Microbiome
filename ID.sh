#!/bin/bash

#SBATCH --time=14-00:00:00   # walltime
#SBATCH -p defq
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH --array=0-3           # job array of size 77
SAVEIFS=$IFS   # Save current IFS
IFS=$'\n'      # Change IFS to new line
myList=('14A'
'14C'
'21D'
'23A')

filename=${myList[${SLURM_ARRAY_TASK_ID}]}
pathoscope ID -alignFile results_"$filename"/"$filename".sam -fileType sam -outDir results_"$filename" -expTag "$filename"
$IFS=$SAVEIFS
