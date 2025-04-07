#!bin/bash

#Dependencies
# StringTie2
# FEELnc
# CPAT
# gffcompare
# bedtools
# Python dependencies - Bio, pandas, seaborn
# HL annotation toolkit

# This pipeline uses FEELnc tool and CPAT to annotate lncRNAs. We will need the following files as input
# 1. BAM files generated from splice aware aligners (Hisat2 or Minimap)
# 2. Reference annotation in .gtf format (make sure all annotation files have the same name across all folders)
# 3. reference genome (make sure there is no other fasta file in the run folder)
# 4. Hexamer and logit files to run CPAT

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <name_of_reference_annotation_file> <name_of_ref_genome_file>"
  exit 1
fi

# Assign the provided arguments to variables
REF_GTF="$1"
REF_GENOME="$2"

# List all of the folders that needs to be processed for the lncRNA pipeline

RUN_FOLDERS=('Bisbetu')


for FOLDER in ${RUN_FOLDERS[@]}; 
do
  cd $FOLDER
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Processing $FOLDER"
  
  # Gather all the bam files into a list
  BAM_FILES=(*.bam)
  
  echo "The following BAM files will be run through StringTie2"
  echo "${BAM_FILES[@]}"
  
  # Running StringTie2 to assemble tanscripts from short reads 
  # Make sure that all of the annotation files (.gtf) files in the folders that need to be processed are having the same name (braker.gtf in this case). 

  OUTPUT_SUFFIX=$FOLDER
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running StringTie"
  stringtie ${BAM_FILES[@]} -G $REF_GTF -o $OUTPUT_SUFFIX.stringtie_transcripts.gtf -p 12 -v 2> $OUTPUT_SUFFIX.stringtie_run.log
  
  # FEElnc run
  # Create a directory for the FEElnc run
  
  # We need to initialise conda in this script
  source "$(conda info --base)/etc/profile.d/conda.sh"
  
  conda activate /home/ae774/feelnc_install_dir

  mkdir FEELnc_run
  cd FEELnc_run
  
  # We need to extract only the transcripts from the reference annotation
  # awk '$3 == "exon" || $3 == "CDS"' ../$REF_GTF > $OUTPUT_SUFFIX.reference_transcripts.gtf
  
  # Perform the FEELnc filter step
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running FEELnc filter"
  FEELnc_filter.pl -i ../$OUTPUT_SUFFIX.stringtie_transcripts.gtf \
  -a ../$REF_GTF \
  -b transcript_biotype="protein_coding" \
  -o $OUTPUT_SUFFIX.stringtie_transcripts.feelncfilter.log > $OUTPUT_SUFFIX.FEELnc_filter_candidate_lncRNA.gtf
  #  -f 0.2 
  
  # Run FEELnc_codpot step
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running FEELnc codpot"
  FEELnc_codpot.pl -i $OUTPUT_SUFFIX.FEELnc_filter_candidate_lncRNA.gtf \
  -a ../$REF_GTF \
  -g ../$REF_GENOME --mode=shuffle 2> $OUTPUT_SUFFIX.feelnc_codpot.log

  conda deactivate
  
  cd ..
  
  # Manual filtering before CPAT run
  
  mkdir CPAT
  cd CPAT
  mkdir manual_filtering
  cd manual_filtering
  
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running Manual filtering step"
  
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Filtering mono-exonic transcripts"
  # STEP 1 - Filtering out the mon0 exonic transcripts
  awk '$3=="exon"' ../../$OUTPUT_SUFFIX.stringtie_transcripts.gtf | \
  awk '{match($0, /transcript_id "([^"]+)"/, a); print a[1]}' | \
  sort | uniq -c | awk '$1 >=2 {print $2}' > more_than_2_exons.txt
  
  grep -F -w -f <(sed 's/^/transcript_id "/; s/$/";/' more_than_2_exons.txt) \
  ../../$OUTPUT_SUFFIX.stringtie_transcripts.gtf > $OUTPUT_SUFFIX.transcripts_exon_filtered.gtf
  
  # STEP 2 - Filtering overlapping transcripts
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Filtering reference overlapping transcripts"
  # Creating bed files with just the exon coordinates, make sure the 12th/10th column has transcript_id in stringtie and reference gtf
  awk '$3 == "exon" {print $1, $4, $5, $12}' OFS="\t" $OUTPUT_SUFFIX.transcripts_exon_filtered.gtf > $OUTPUT_SUFFIX.candidate_exons.bed
  awk '$3 == "exon" {print $1, $4, $5, $14}' OFS="\t" ../../$REF_GTF > $OUTPUT_SUFFIX.reference_exons.bed # USER INPUT required!
  
  bedtools intersect -v -a $OUTPUT_SUFFIX.candidate_exons.bed -b $OUTPUT_SUFFIX.reference_exons.bed | cut -f4 > transcripts_no_overlap.txt
  
  awk '{match($0, /"([^"]+)"/, a); print a[1]}' transcripts_no_overlap.txt| sort | uniq > temp.txt
  rm transcripts_no_overlap.txt
  mv temp.txt transcripts_no_overlap.txt
  
  grep -F -w -f <(sed 's/^/transcript_id "/; s/$/";/' transcripts_no_overlap.txt) \
  $OUTPUT_SUFFIX.transcripts_exon_filtered.gtf > $OUTPUT_SUFFIX.manual_filtered_candidate_lncrna.gtf
  
  # Extracting the transcripts into a fasta file
  python /home/ae774/Projects/Tools/annotation_toolkit/run.py -f transcript \
  -r ../../$REF_GENOME $OUTPUT_SUFFIX.manual_filtered_candidate_lncrna.gtf
  
  cp Annotation_toolkit_outputs/transcript.fa ../$OUTPUT_SUFFIX.manual_filtered_candidate_lncrna.fa
  
  cd ..
  
  # CPAT run
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running CPAT"
  # We will need the hexamer and the logit files pre-bilt before running the pipeline
  HEXAMER_FILE_PATH="/home/ae774/Projects/lncRNA/moth_melanism/RNA_seq/Bisbetu/bisbetu_hexamer.tsv" # USER INPUT
  LOGIT_FILE_PATH="/home/ae774/Projects/lncRNA/moth_melanism/RNA_seq/Bisbetu/bisbetu.logit.RData" # USER INPUT
  
  cpat -x $HEXAMER_FILE_PATH -d $LOGIT_FILE_PATH --top-orf=100  --antisense \
  -g $OUTPUT_SUFFIX.manual_filtered_candidate_lncrna.fa \
  -o $OUTPUT_SUFFIX.cpat_output
  
  # we need to combine the transcripts with no ORFs and the ones with really low coding probability
  
  # Extracting the transcripts with coding probability < 0.4
  awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Coding_prob") cp=i; else if ($i == "ID") id=i; next} $cp < 0.4 {print $id}' \
  $OUTPUT_SUFFIX.cpat_output.ORF_prob.tsv > coding_prob_0.4.txt
  
  sed 's/_ORF_[0-9]\+//' coding_prob_0.4.txt | sort -u > coding_prob_0.4_clean.txt
  
  cat coding_prob_0.4_clean.txt $OUTPUT_SUFFIX.cpat_output.no_ORF.txt > final_lncRNA_transcripts.txt
  
  grep -F -w -f <(sed 's/^/transcript_id "/; s/$/";/' final_lncRNA_transcripts.txt) \
  manual_filtering/$OUTPUT_SUFFIX.manual_filtered_candidate_lncrna.gtf > $OUTPUT_SUFFIX.CPAT_final_lncRNA.gtf
  
  cd ..

  # Merging the annotations
  
  # creating a folder and copying the final sets of lncRNAs from the FEELnc and the CPAT methods
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Merging the FEELnc and the CPAT results"
  mkdir merger
  cd merger
  cp ../FEELnc_run/feelnc_codpot_out/$OUTPUT_SUFFIX.FEELnc_filter_candidate_lncRNA.gtf.lncRNA.gtf .
  cp ../CPAT/$OUTPUT_SUFFIX.CPAT_final_lncRNA.gtf .
  
  gffcompare $OUTPUT_SUFFIX.FEELnc_filter_candidate_lncRNA.gtf.lncRNA.gtf $OUTPUT_SUFFIX.CPAT_final_lncRNA.gtf
  
  mv gffcmp.combined.gtf $OUTPUT_SUFFIX.merged_final_consensus_lncRNAs.gtf
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Final set of lncRNAs are in  $OUTPUT_SUFFIX.merged_final_consensus_lncRNAs.gtf"
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Finished processing $OUTPUT_SUFFIX"
  echo "----------------------------DONE----------------------------------"
  cd ..
  
  cd ..
done




















