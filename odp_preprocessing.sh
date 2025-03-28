# !bin/bash
# This script is used to clean the set of lncRNA annotations obtained from the HL lncRNA pipeline

# We need the following files
# 1. fasta file of annotated lncRNAs
# 2. gtf file 

# Following dependencies are needed
# seqkit
# EMboss transeq


# Declare the input samples
samples=('011-UU' '053-LA' '062-UU' '070-UU' '081-UU' '082-UU')

for sample in ${samples[@]};
do
  # converting the fasta into .pep file
  transeq -sequence $sample.merged_final_consensus_lncRNAs.fa -outseq temp.pep -clean
  
  # transeq adds in an _ to each transcript_id tha we need to remove using the below command
  awk '/^>/ {sub("_1$", "");} {print}' temp.pep > temp1.pep

  # Identifying and removing the duplicate entried from the fasta file
  seqkit rmdup -s -i -o $sample.merged_final_consensus_lncRNAs_nodups.pep \
  temp1.pep -d $sample.duplicates.fa
  
  # Take the duplicate transcript IDs into a separate file
  grep '^>' $sample.duplicates.fa | sed 's/^>//' > $sample.duplicates.txt
  
  # create the chrom file from the .gtf file
  awk '$3 == "transcript" {print $10, $1, $7, $4, $5}' \
  ../../Samples/$sample/merger/$sample.merged_final_consensus_lncRNAs.gtf | uniq > temp1.chrom
  
  # clean the temp1.chrom by removing " and ; from from the ID column
  awk '{gsub(/[";]/, "", $1); print}' temp1.chrom > temp2.chrom
  
  # The above only gives a space separated file and we need to convert that into a tab separated file
  tr -s ' ' '\t' < temp2.chrom > temp3.chrom
  
  # removing the duplicates from the chrom file
  grep -v -F -f $sample.duplicates.txt temp3.chrom > $sample.merged_final_consensus_lncRNAs.chrom
  
  # removing the duplicates from the .gtf
  grep -v -F -f <(sed 's/^/transcript_id "/; s/$/";/' $sample.duplicates.txt) \
  ../../Samples/$sample/merger/$sample.merged_final_consensus_lncRNAs.gtf > $sample.merged_final_consensus_lncRNAs_nodups.gtf
  
  # removing the temp files
  rm temp*
  
done
 
  
  