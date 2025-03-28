# !bin/bash

# Run featue counts for all of the samples

samples=('011-UU' '053-LA' '062-UU' '070-UU' '081-UU' '082-UU')

for sample in ${samples[@]};
do
  
  echo "Counting $sample"
  bam_files=(../$sample/*.bam)
  echo "Features counted using the following bam files ${bam_files[@]}"
  
  featureCounts -a ../$sample.merged_final_consensus_lncRNAs_nodups.gtf \
  -o $sample.feature_counts.tsv ${bam_files[@]} \
  -F 'GTF' -t transcript -g transcript_id \
  -M -T 12 -f -p
  
done