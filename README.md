# undergraduate_level_cfDNA_lab
Generate simulated "sam-like" files for teaching students how to handle cell-free genomic data

# Intention
This lab was created to help students familiarize themselves with genomic data. The scripts here create "Sam-like" files that
are tsv's with READID , CHR, START, END, TLEN, although any sam column can be added.

# How to use
The script will take in a genome file and index it if needed. This genome file can be something like hg19, or hg19 concatenated
with  microbial genome if infection is something you want to study.

# Diseases simulated
1- Aneuploidy in pregnant women
2- Transplant rejection
3- Bloodstream infection

The data generated has length and fragmentation patterns typically found in cfDNA experiments

# Use
```
python scripts/generate_sam.py --num_reads 1000000 \
                              --genome_file PATH_TO_FASTA \
                              --length_distribution_file PATH_TO_FILE \
                              --patient_sex (M or F) \
                              --output_fle (outfile) \
                              --BSI_reads (number of reads for bloodstrem infection simulation) \
                              --donor_sex (M or F) \
                              --donor_frac (float) \
                              --chr_21_fraction (float) \
```
