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

#/home/users/nus/e1101919/.conda/envs/muscle/bin/Rscript /home/users/nus/e1101919/scratch/software/leafcutter/scripts/leafcutter_ds.R \
#							--num_threads 23 -t 300 -o diff_leafcutter/leafcutter_ds.sex.${c} aida_perind_numers.${c}.counts.gz groups_file.sex.${c}.txt

for n in `seq 8`
do
    /home/users/nus/e1101919/.conda/envs/muscle/bin/Rscript /home/users/nus/e1101919/scratch/software/leafcutter/scripts/leafcutter_ds.R \
							    --num_threads 23 -t 300 -o diff_leafcutter/leafcutter_ds.sex.${n}PC.${c} aida_perind_numers.${c}.counts.gz groups_files/groups_file.sex.${n}PC.${c}.txt
done
