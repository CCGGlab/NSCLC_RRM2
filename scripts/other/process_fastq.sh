#!/bin/bash

# Get all files to be processed
selected_files=$(ls raw/NSCLCL/last_sample/X204SC21041081_Z01_F002/rawdata/*1.fq.gz)

# Process one by one
for f in $selected_files;
do
        echo 'Processing' $f;
  f="${f##*/}"
  file1=$f
  file2="${file1/1.fq.gz/2.fq.gz}"

  # Check whether processed already and skip if it does
  if [ -f raw/gene_counts/$f'.gene_counts' ];
  then
      echo "raw/gene_counts/$f'.gene_counts' exists already, skipping ...";
      continue;
  fi;

  # Align to grch38 using hisat: +/-1h
        echo '  Aligning ...';
  (hisat2 -p 8 --no-mixed --no-discordant --rna-strandness RF --no-unal -x downloads/genomes/igenomes/hg38/genome -1 raw/NSCLCL/last_sample/X204SC21041081_Z01_F002/rawdata/$file1 -2 raw/NSCLCL/last_sample/X204SC21041081_Z01_F002/rawdata/$file2 | samtools view -bS - > raw/$f'.bam') >& raw/hisat2_log/$f'_log.txt'

  # Quantify to gencode 29 using htseq-count: +/-2h30
        echo '  Quantifying ...';
  (samtools view raw/$f'.bam' | htseq-count -m intersection-strict -s reverse -i gene_id - downloads/gencode/gencode.v29.annotation.gtf > raw/gene_counts/$f'.gene_counts') >& /dev/stdout | tail -n200 > raw/htseq_log/$f'_log.txt'
  # samtools view -H raw/AS1h_1.bam # "chr1" format, "1" format for igenomes grch38

  # Get bam stats: includes insert sizes, read lengths, ...
  samtools stats raw/$f'.bam' > raw/bam_stat/$f'_stat_log.txt'

  # Rm bam file
        echo '  Finishing ...';
  rm raw/$f'.bam'
  echo ''

done

