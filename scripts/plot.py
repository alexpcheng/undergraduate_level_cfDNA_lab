#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

def parse_args(chromosomes):
    file = sys.argv[1]
    df = pd.read_csv(file, sep = '\t')

    genome_sizes = {}
    with open(sys.argv[2]) as f:
        for line in f:
            key, val = line.strip().split('\t')[0:2]
            if key in chromosomes:
                genome_sizes[key] = int(val)

    outfile = file.replace('.txt', '.pdf')
    return(df, genome_sizes, outfile)


def get_fractions(df, chr_lengths, chromosomes):
    total_genome_size=0
    for key, val in chr_lengths.items():
        total_genome_size+=val

    sequenced_reads = {}
    effective_reads = {}
    total_reads = 0
    effective_total_reads = 0
    for chr in chromosomes:
        reads = len(df[df['CHR']== chr].index)
        eff_reads = reads / (chr_lengths[chr] / total_genome_size)
        total_reads += reads
        effective_total_reads+= eff_reads
        sequenced_reads[chr]=reads
        effective_reads[chr]=eff_reads

    sequenced_fraction = {}
    effective_fraction = {}
    for chr in chromosomes:
        seq_frac = sequenced_reads[chr]/total_reads
        eff_frac = effective_reads[chr]/effective_total_reads
        sequenced_fraction[chr] = seq_frac
        effective_fraction[chr] = eff_frac

    chr1_fraction = effective_fraction['chr1']
    chrX_fraction = effective_fraction['chrX']
    chrY_fraction = effective_fraction['chrY']
    BSI_fraction = effective_fraction['NC_000913.2']
    chr21_fraction = effective_fraction['chr21']
    male_frac = 2*effective_fraction['chrY']/effective_fraction['chr1']
    female_frac = 1-male_frac

    return(chr1_fraction, chrX_fraction, chrY_fraction, BSI_fraction, chr21_fraction, male_frac, female_frac)


def plot(df, chr1_fraction, chrX_fraction, chrY_fraction, BSI_fraction, chr21_fraction, male_frac, female_frac, outfile):

    lengths = pd.to_numeric(df['TLEN'], errors= 'coerce')
    print(lengths)
    #print(lengths)
    lengths = lengths[lengths >=0]
    d = pd.DataFrame(lengths.value_counts(normalize=True, sort=True))
    d['length'] = d.index
    d = d.sort_values(by=['length'], axis=0)
    plt.plot(d['length'], d['TLEN']*100)
    ann = f'chr1: {chr1_fraction:.5f}\nchr21: {chr21_fraction:.5f}\nchrX: {chrX_fraction:.5f}\nchrY: {chrY_fraction:.5f}\nBSI fraction: {BSI_fraction:.5f}\nFemale fraction: {female_frac:.5f}\nMale fraction: {male_frac:.5f}'
    plt.text(s = ann, x=325, y=1.25)
    plt.xlabel('Fragment length (base pairs)')
    plt.ylabel('Fraction of reads (%)')
    plt.savefig(outfile)


def main():
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'NC_000913.2']

    df, chr_lengths, outfile = parse_args(chromosomes)
    chr1_fraction, chrX_fraction, chrY_fraction, BSI_fraction, chr21_fraction, male_frac, female_frac = get_fractions(df, chr_lengths, chromosomes)
    plot(df, chr1_fraction, chrX_fraction, chrY_fraction, BSI_fraction, chr21_fraction, male_frac, female_frac, outfile)

if __name__ == '__main__':
    main()
