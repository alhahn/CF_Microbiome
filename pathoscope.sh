#!/bin/bash

#SBATCH --time=14-00:00:00   # walltime
#SBATCH -p defq
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH --array=0-3           # job array of size 98

SAVEIFS=$IFS   # Save current IFS
IFS=$'\n'      # Change IFS to new line
 
myList=('14A'
'14C'
'23A'
'21D')
filename=${myList[${SLURM_ARRAY_TASK_ID}]}
pathoscope MAP -1 ../"$filename"_R1_001.fastq -2 ../"$filename"_R2_001.fastq -numThreads 20 -expTag "$filename" -indexDir ./index_"$filename" -targetRefFiles AllBacGenomes.fna -outDir ./results_"$filename" -outAlign "$filename".sam
$IFS=$SAVEIFS
