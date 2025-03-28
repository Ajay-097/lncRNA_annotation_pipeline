# lncRNA_annotation_pipeline
Annotation of lncRNAs from RNA seq data using complementary approaches
# How to run lncRNA pipeline?
The pipeline is designed to annotate lncRNAs by iderating through the sample folders.
The script file needs to be opened and any sample that needs to be processed through the pipeline should be added to the RUN_FOLDERS array list.
Make sure the following files are present in the sample folder
  1. BAM files 
  2. reference annotation in .gtf format (if running multiole folders make sure the annotation file name is same across all the folders)
  3. Reference genome with .fasta suffix (make sure there is no other . fasta files)
  4. Hexamer and logit files to run CPAT
Run the script from the directory where all smaple folders are present  
