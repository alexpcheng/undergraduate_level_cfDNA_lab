#!/usr/bin/env bash
tail -n +2 synthetic_samples.tsv | while read line
do
  outfile=$(echo $line | cut -d' ' -f1)
  pt_sex=$(echo $line | cut -d' ' -f2)
  donor_sex=$(echo $line | cut -d' ' -f3)
  BSI_frac=$(echo $line | cut -d' ' -f4)
  chr21_frac=$(echo $line | cut -d' ' -f5)
  donor_frac=$(echo $line | cut -d' ' -f6)
  num_reads=$(echo $line | cut -d' ' -f7)

  if [[ -f $outfile ]]; then
    echo "$outfile exists. will not overwrite"
  else
    echo ">>> Creating $outfile <<<"

    python scripts/generate_sam.py \
                --num_reads $num_reads \
                --genome_file references/genome.fa \
                --length_distribution_file references/length_dist.txt \
                --patient_sex $pt_sex \
                --donor_sex $donor_sex \
                --BSI_reads $BSI_frac \
                --chr21_fraction $chr21_frac \
                --donor_frac $donor_frac \
                --output_file $outfile
    python scripts/plot.py $outfile references/genome.fa.fai
  fi
done
