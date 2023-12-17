#!/bin/bash

# c6i.12xlarge

# running in screen (avoid being interupted by broken shell)
#screen -S psb
#./jobs_

mkdir -p /ebs2/aida-sqtl/FreezeV1_0/03_pseudobulk_bam/meta
cat ../01_STARsolo/list_sample.txt | while read line
do
    # get names firstly
    fullname=`echo $line | cut -d' ' -f 1`
    abbrname=`echo $line | cut -d' ' -f 2`

    batch=`echo $fullname | sed 's/.*B0*\([0-9]*\)_L00.*/\1/'`
    library=`echo $fullname | sed 's/.*L00*\([1-2]*\)_5GEX.*/\1/'`
    country=`echo $fullname | cut -c 1-2`

    lib_name=`echo ${country}_B${batch}_L${library}`


    # prepare processed file (bam + meta) 
    echo "pseudobulk for $fullname"
    ls -l /ebs2/aida-sqtl/FreezeV1_0/02_dedup/${fullname}_dedup.bam

#    if [ ! -f /ebs2/aida-sqtl/FreezeV1_0/03_pseudobulk_bam/meta/meta_$fullname.txt ]
#    then
    # col9 == $lib_name & col11 == $abbr_name
    less /ebs1/users/tongyihan/AIDA/data/14.individual_cell_num/AIDA_DataFreeze_v1_2_cell_annotations_three_levels_general_flags_merge.txt |\
	grep $lib_name | grep $abbrname | cut -f1,3 | sed "s/-$lib_name//" > /ebs2/aida-sqtl/FreezeV1_0/03_pseudobulk_bam/meta/meta_$fullname.txt
#    fi


    # start processing (pseudobulk)
    nohup python pseudobulk_bam.py \
	  --bam_file /ebs2/aida-sqtl/FreezeV1_0/02_dedup/${fullname}_dedup.bam \
	  --csv_file /ebs2/aida-sqtl/FreezeV1_0/03_pseudobulk_bam/meta/meta_$fullname.txt \
	  --out_prefix /ebs2/aida-sqtl/FreezeV1_0/03_pseudobulk_bam/${fullname} > /s3/aida-sqtl/FreezeV1_0/00_run_log/04_pseudobulk_$fullname.log 2>&1 &


    # avoid submit too quickly (pgrep may delay), could be set to 1000 when testing using the first sample
    sleep 1
    

    # check for number of running jobs constantly, keep running # jobs
    while true
    do
	numJobs=`pgrep -cx python -U tianchi`
#	echo -e "checking # jobs: $numJobs"
	if [ $numJobs -lt 15 ]
	then
	    break
	fi
	
	sleep 1
    done
done
