import time
from tqdm import tqdm
import argparse
import pysam
import os
import random

from numba import jit

celltype = dict()
seq_list = dict()

def read_bam(bam_file, out_prefix):
    t0 =  time.time()
    print("Start to process...")
    xbarcode_reads = 0
    multimap_reads = 0
    waspfilt_reads = 0
    retained_reads = 0
    wasp_reads = 0
    n_print = 0 

    # load bam file
    input_file = pysam.AlignmentFile(bam_file,'rb')

    output_files = {}
    for ct in list(set(celltype.values())):
        if os.path.exists('.'.join( (out_prefix,ct,'bam') )):
            os.remove('.'.join( (out_prefix,ct,'bam') ))
        output_files[ct] = pysam.AlignmentFile( '.'.join( (out_prefix,ct,'bam') ),'wb', template=input_file)
        print('.'.join( (out_prefix,ct,'bam')))
#    all_reads = input_file.read().split('\n')
#    len_reads = len(all_reads)
#    t1 = time.time()
#    print("load bam file successfully, time consumption:", t1-t0, "reads number is:", int(len_reads))

#    for i in tqdm(input_file):
    for i in input_file:
        barcode = i.get_tag('CB')

        # barcode annotated 
        if barcode in celltype.keys():

            # 1. UMI collapse not applicable (R2 of same UMI can cover different exons)
            # Instead, Picard could be used to remove duplicates, which consist of both the same R1 and the sanme R2
            
            # 2. wasp filter
            try:
                wasp = i.get_tag('vW')                    
            except:
                Nmap = i.get_tag('NH')
                if Nmap == 1:
                    retained_reads += 1
                    output_files[celltype[barcode]].write(i)
                else:
                    multimap_reads += 1
            else:
                wasp_reads += 1
                if wasp == 1:
                    retained_reads += 1
                    output_files[celltype[barcode]].write(i)
                else:
                    waspfilt_reads += 1
        else:
            xbarcode_reads += 1
                
    input_file.close()
    for ct in list(set(celltype.values())):
        output_files[ct].close()

    t2 = time.time()
    print("convert sucessfully, total time:", t2-t0)

    print("xbarcode_reads:", xbarcode_reads)
    print("multimap_reads:", multimap_reads)
    print("waspfilt_reads:", waspfilt_reads)
    print("retained_reads:", retained_reads)
    print("wasp_reads:", wasp_reads)


def read_csv(csv_file):
    with open(csv_file) as f:
        for line in f:
            line_dict = line.split('\n')[0].split('\t')
            celltype[line_dict[0]] = line_dict[1]

def main():
    parser = argparse.ArgumentParser(description='Set file path')
    parser.add_argument('--bam_file', help='source bam file')
    parser.add_argument('--csv_file', help='source meta file')
    parser.add_argument('--out_prefix', help='outfile prefix')

    args = parser.parse_args()

    read_csv(args.csv_file)
    read_bam(args.bam_file, args.out_prefix)

    

if __name__ == "__main__":

    main()
