#!/bin/bash


grep -l "finished successfully" /ebs2/aida-sqtl/FreezeV1_0/01_STARsolo/*/log_files/*.out | sed 's/.*log_files\/\(.*\)\.out/\1/' | grep -v -f - list_sample.6-21.txt > list_sample.6-21.remain.txt 
## run by

# 1. keep running in background (nohup)
#bash batch_run.1by1.sh

# 2. running in screen (avoid being interupted by broken shell)
#screen -S star
#./jobs_02_STARsolo.1by1.sh 


cat list_sample.6-21.remain.txt | while read sample
do
    echo "STAR mapping for $sample"
    fullname=`echo $sample | cut -d' ' -f1`
    abbrname=`echo $sample | cut -d' ' -f2`
    echo $fullname, $abbrname
    nohup bash 02_run_STARsolo.sh $fullname $abbrname > /s3/aida-sqtl/FreezeV1_0/00_run_log/02_STARsolo_$fullname.log 2>&1 &

    # avoid submit too quickly (pgrep may delay)
    sleep 3

    # check for number of running jobs constantly, keep running 15 jobs
    # 500GB RAM, 64 vCPU: 15 jobs (each ~ 30GB RAM, 4 threads)
    while true
    do
	numJobs=`pgrep -cx STAR`
#	echo -e "checking # jobs: $numJobs"
	if [ $numJobs -lt 28 ]
	then
	    break
	fi
	
	sleep 3
    done
done
