#!/bin/bash

# Shell to use
#$ -S /bin/bash
# Name of the job in SGE
#$ -N my_job
# Name of the queue to use
#$-q cemeb.q
# Maximum hardware time allowed for this job
#$-l h_rt=15:00:00
# Merge output (error and standard with '-j y') into one single file called my_job.out
#$-o my_job.out -j y
# run in the current directory
#$ -cwd

cd /home/rvezy/2-Simulations_FullPeriod2012_Seq
echo "$SGE_TASK_ID" | ./maespa
tar czf  ${SGE_TASK_ID}_maespa.tar.gz ${SGE_TASK_ID}_hrflux.dat ${SGE_TASK_ID}_Dayflx.dat ${SGE_TASK_ID}_Maeserr.dat ${SGE_TASK_ID}_watbal.dat ${SGE_TASK_ID}_watbalday.dat ${SGE_TASK_ID}_watlay.dat ${SGE_TASK_ID}_watsoilt.dat ${SGE_TASK_ID}_watupt.dat
rm ${SGE_TASK_ID}_hrflux.dat
rm ${SGE_TASK_ID}_Dayflx.dat
rm ${SGE_TASK_ID}_Maeserr.dat
rm ${SGE_TASK_ID}_watbal.dat
rm ${SGE_TASK_ID}_watbalday.dat
rm ${SGE_TASK_ID}_watlay.dat
rm ${SGE_TASK_ID}_watsoilt.dat
rm ${SGE_TASK_ID}_watupt.dat
