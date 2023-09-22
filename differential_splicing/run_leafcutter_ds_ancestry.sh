#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=23:mem=96G
#PBS -l walltime=24:00:00
#PBS -P 11003054
#PBS -N leafcutter_ds
#PBS -o leafcutter_ds.$c.o
#PBS -e leafcutter_ds.$c.e



cd $PBS_O_WORKDIR

mkdir -p diff_leafcutter
ls group_file/ | grep $c.txt | while read f
do
    v=`echo $f | sed 's/.*\.country\.\(.*\)\..*/\1/'`
    #echo $v
    /home/users/nus/e1101919/.conda/envs/muscle/bin/Rscript /home/users/nus/e1101919/scratch/software/leafcutter/scripts/leafcutter_ds.R \
							    --num_threads 23 -t 300 -o diff_leafcutter/leafcutter_ds.country.${v} ../DSG_leafcutter_sex/aida_perind_numers.${c}.counts.gz group_file/groups_file.country.${v}.txt
done
