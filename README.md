# scRNASeq_ZSWIM8_public_dataset_MK
location:/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al
location:/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/CO_DiO_DAY100_ScRNA/20250212_AV241602_2_12_2025_A_Manoj/Samples/Manoj

# 150 *2 bp
#aviti element
#20,000 cells 

# md5sum check for integrity 
md5sum -c 20250212_AV241602_2_12_2025_A_Manoj.md5 
CO_DAY100_R1.fastq.gz: OK
CO_DAY100_R2.fastq.gz: OK

#https://nf-co.re/scrnaseq/3.0.0/
 eval "$(/home/sor4003/anaconda3/bin/conda shell.bash hook)"
 conda activate env_nf
 nextflow pull nf-core/scrnaseq -r 3.0.0
