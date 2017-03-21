#!/bin/bash

# Shell to use
#$ -S /bin/bash
# Name of the job in SGE
#$ -N Euca_2012
# Name of the queue to use
#$-q cemeb.q
# Maximum hardware time allowed for this job
#$-l h_rt=30:00:00
# Merge output (error and standard with '-j y') into one single file called my_job.out
#$-o my_job.out -j y
# run in the current directory
#$ -cwd

echo "$SGE_TASK_ID" | maespa
