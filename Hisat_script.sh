#!/bin/bash
# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <path_to_samples> <path_to_genome_index>"
  exit 1
fi

# Assign the provided arguments to variables
SAMPLES_PATH="$1"
HISAT_INDEX="$2"

for R1 in "$SAMPLES_PATH"/*_1.fastq;
do
  # Run hisat for each of the read files
  R2=${R1//_1.fastq/_2.fastq}
  OUTPUT_FILE=${R1//_1.fastq/output}
  OUTPUT_FILE=${OUTPUT_FILE##*/}
  echo "R1 - $R1"
  echo "R2 - $R2"
  echo "Output is stored to $OUTPUT_FILE"
  hisat2 -x $HISAT_INDEX -1 $R1 -2 $R2 --dta -p 12 -S $OUTPUT_FILE.sam 2> $OUTPUT_FILE.hisat.log
  
  # sorting and conversion to bam file
  samtools sort -@ 12 $OUTPUT_FILE.sam > $OUTPUT_FILE.sorted.bam  
  
  # remove the sam file and un sorted bam file
  # rm $OUTPUT_FILE.sam
done