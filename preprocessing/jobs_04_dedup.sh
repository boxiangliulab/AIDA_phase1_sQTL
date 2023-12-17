#!/bin/bash

# m6i.24xlarge	

# running in screen (avoid being interupted by broken shell)
#screen -S psb
#./jobs_

mkdir -p /ebs2/aida-sqtl/FreezeV1_0/02_dedup

# should use all samples in ../01_STARsolo/list_sample.txt 
#grep -f ../summary_sample_use.txt ../01_STARsolo/list_sample.txt | while read line
#grep -v -f ../summary_sample_use.txt ../01_STARsolo/list_sample.txt | while read line
cat ../01_STARsolo/list_sample.txt | while read line
do
    # get names firstly
    fullname=`echo $line | cut -d' ' -f 1`
    abbrname=`echo $line | cut -d' ' -f 2`


    echo "deduplicate for $fullname"
    ls -l /s3/aida-sqtl/FreezeV1/01_STARsolo/$abbrname/${fullname}_Aligned.sortedByCoord.out.bam
    # start processing (deduplicate)
    # Exception in thread "main" htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found: A00552:133:H5KTTDSX2:4:1311:31439:25191 2/2 151b aligned to chr1:10002-10113.
    # We care about both mapped position and UMI
    nohup java -jar /ebs1/shared/software/picard.jar MarkDuplicates \
	  REMOVE_DUPLICATES=true \
	  BARCODE_TAG=CB \
	  MOLECULAR_IDENTIFIER_TAG=UB \
	  I=/s3/aida-sqtl/FreezeV1/01_STARsolo/$abbrname/${fullname}_Aligned.sortedByCoord.out.bam \
	  O=/ebs2/aida-sqtl/FreezeV1_0/02_dedup/${fullname}_dedup.bam \
	  M=/ebs2/aida-sqtl/FreezeV1_0/02_dedup/marked_dup_metrics.$fullname.txt > /s3/aida-sqtl/FreezeV1_0/00_run_log/03_dedup_$fullname.log 2>&1 &
    
    
    # avoid submit too quickly (pgrep may delay), could be set to 1000 when testing using the first sample
    sleep 1

    # check for number of running jobs constantly, keep running # jobs
    while true
    do
	numJobs=`pgrep -cx java -U tianchi`
	# https://github.com/hzi-bifo/phylogeny-of-single-cells/commit/a3601bf029fb0545984da667dcfb13d0388bee38
	# 30G x 
	if [ $numJobs -lt 10 ]
	then
	    break
	fi
	
	sleep 1
    done
done
