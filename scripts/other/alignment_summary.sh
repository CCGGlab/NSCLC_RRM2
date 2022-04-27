#!/bin/bash

# Check alignment efficiency
############################

# Create new  file with summary
echo -n > results/data/alignment_summary.txt

# Extract all efficiencies
selected_files=$(ls -d -1 $PWD/raw/hisat2_log/*)

# Process one by one
for f in $selected_files; 
do 
# echo 'Processing' $f;
echo -n -e $f "\t" >> results/data/alignment_summary.txt;
egrep -i '\<alignment rate\>' $f | grep -o '[0-9]*\.[0-9]*' >> results/data/alignment_summary.txt
done
