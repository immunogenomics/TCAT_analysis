#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.tsv>"
    exit 1
fi

# Assign command line arguments to variables
input_file=$1
output_file=$2


samtools view -q 255 "$input_file" | \
awk 'BEGIN {OFS="\t"} 
     $0 ~ /CB:/ && $0 ~ /RE:/ && $0 ~ /UB:/ {
         for(i=1;i<=NF;i++) {
             if($i ~ /CB:/) cb=$i;
             else if($i ~ /RE:/) re=$i;
             else if($i ~ /UB:/) ub=$i;
         }
         print cb, re, ub
     }' > "$output_file"