#!/usr/bin/env python3

import sys
from pyfaidx import Fasta
import argparse
import random
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--num_reads", help = "The approximate number of reads to generate")
    parser.add_argument("--genome_file", help = "The genome file to use")
    parser.add_argument("--length_distribution_file", help = "File length distribution file")
    parser.add_argument("--patient_sex", help = "Sex of the person (for chromosome X and Y)")
    parser.add_argument("--output_file")
    parser.add_argument("--BSI_reads", type = int, help = "Bloodstream infection reads", default = 0)
    parser.add_argument("--donor_sex", help = "Sex of donor, implies transplant patient", choices = ['male', 'female'])
    parser.add_argument("--donor_frac", type= float,  default = 0)
    parser.add_argument("--chr21_fraction", type= float, help = "Chr21 fraction, implies pregnant patient")

    args = parser.parse_args()

    num_reads = int(args.num_reads)
    genome_file = args.genome_file
    len_file = args.length_distribution_file
    patient_sex = (args.patient_sex).lower()
    outfile = args.output_file

    BSI_reads = args.BSI_reads
    donor_sex = (args.donor_sex).lower()
    donor_frac = float(args.donor_frac)
    chr21_frac = args.chr21_fraction

    l1 = []
    l2 = []
    with open(len_file) as f:
        for line in f:
            a, b = line.strip().split('\t')
            l1.append(float(a))
            l2.append(int(b))
    lengths = [l2, l1]
    return(num_reads, genome_file, lengths, patient_sex, BSI_reads, donor_sex, donor_frac, chr21_frac, outfile)

def load_reference_genome(genome_file, chromosomes):
    genome = Fasta(genome_file) #generates index if it doesnt exist
    genome_sizes = {}
    with open(genome_file + '.fai') as f:
        for line in f:
            key, val = line.strip().split('\t')[0:2]
            if key in chromosomes:
                genome_sizes[key] = int(val)
    return(genome, genome_sizes)

def reads_per_chromosome(num_reads, genome, genome_sizes, BSI_reads, patient_sex, donor_sex, donor_frac, chr21_frac, chromosomes):
    total_genome_size = sum(genome_sizes.values())

    reads_per_chr = {}
    for chr in chromosomes:
        reads_to_take = round(num_reads/total_genome_size*genome_sizes[chr])
        if chr == 'chrX' or chr == 'chrY':
            reads_to_take = round(reads_to_take/2)
        reads_per_chr[chr] = reads_to_take

    # Adjust based on pathology
    if patient_sex == 'female' and donor_sex != 'male':
        reads_per_chr['chrY'] = 0
        reads_per_chr['chrX'] = 2*reads_per_chr['chrX']

    if (patient_sex == 'male' and donor_sex == 'female'): #would have less Y now ...
        reads_per_chr['chrY'] = round(reads_per_chr['chrY'] - (donor_frac)*reads_per_chr['chrY'])


    if (patient_sex == 'female' and donor_sex == 'male'):
        reads_per_chr['chrY'] = round((donor_frac)*reads_per_chr['chrY'])
        #reads_per_chr['chrX'] = round(reads_per_chr['chrX'] + donor_frac*reads_per_chr['chrX'])


    if (chr21_frac is not None):
        reads_per_chr['chr21'] = round((1+chr21_frac)*reads_per_chr['chr21'])

    if BSI_reads > 0:
        reads_per_chr['NC_000913.2'] = BSI_reads

    return(reads_per_chr)

def get_read(genome_sizes, genome, lengths, chr):
    GOT_Ns = True
    while (GOT_Ns):
        START = random.randint(0, genome_sizes[chr])
        TLEN = np.random.choice(lengths[0], p = lengths[1])
        END = TLEN + START -1
        MOLECULE = str(genome[chr][START:END]).upper()
        if 'N' not in MOLECULE:
            GOT_Ns = False
    R1_strand = np.random.choice(['neg', 'pos'], p = [0.5, 0.5])
    if R1_strand == 'pos':
        START_1 = START
        START_2 = END
        TLEN_1 = TLEN
        TLEN_2 = TLEN
        SEQ_1 = MOLECULE[:75]
        SEQ_2 = MOLECULE[-75:]
    else:
        START_2 = START
        START_1 = END
        TLEN_1 = TLEN
        TLEN_2 = TLEN
        SEQ_2 = MOLECULE[:75]
        SEQ_1 = MOLECULE[-75:]


    return(START_1, START_2, TLEN_1, TLEN_2, SEQ_1, SEQ_2)

def generate_pseudosam(outfile, genome_sizes, genome, reads_per_chr, lengths):
    with open(outfile, 'w') as w:
        count = 1
        w.write('\t'.join(['READ_ID', 'CHR', 'START', 'TLEN', 'SEQ'])+'\n')
        for chr in reads_per_chr:
            for _ in range(0, reads_per_chr[chr]):
                READ_ID_1 = f'read_id_{count:06}_1'
                READ_ID_2 = f'read_id_{count:06}_2'
                CHR = chr
                START_1, START_2, TLEN_1, TLEN_2, SEQ_1, SEQ_2 = get_read(genome_sizes, genome, lengths, chr)
                if chr == 'NC_000913.2':
                    CHR = 'unknown'
                    START_1 = '*'
                    START_2 = '*'
                    TLEN_1 = '*'
                    TLEN_2 = '*'
                string_1 = '\t'.join([READ_ID_1, CHR, str(START_1), str(TLEN_1), SEQ_1])+'\n'
                string_2 = '\t'.join([READ_ID_2, CHR, str(START_2), str(TLEN_2), SEQ_2])+'\n'
                w.write(string_1+string_2)
                count+=1


def main():
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'NC_000913.2']
    num_reads, genome_file, lengths, patient_sex, BSI_reads, donor_sex, donor_frac, chr21_frac, outfile = parse_args()
    genome, genome_sizes = load_reference_genome(genome_file, chromosomes)

    if BSI_reads == 0:
        chromosomes.remove('NC_000913.2')
        del genome_sizes['NC_000913.2']

    reads_per_chr = reads_per_chromosome(num_reads, genome, genome_sizes, BSI_reads, patient_sex, donor_sex, donor_frac, chr21_frac, chromosomes)
    generate_pseudosam(outfile, genome_sizes, genome, reads_per_chr, lengths)
if __name__ == '__main__':
    main()
