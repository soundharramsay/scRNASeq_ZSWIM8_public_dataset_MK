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

# fatsq_naming 
mv CO_DAY135_R1.fastq.gz CO_DAY135_S1_L001_R1_001.fastq.gz
mv CO_DAY135_R2.fastq.gz CO_DAY135_S1_L001_R2_001.fastq.gz


#location of result
/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/1_CO_DiO_DAY100_ScRNA/20250212_AV241602_2_12_2025_A_Manoj/Samples/Manoj/mapping_results_with_fastq

## slurn submission
#!/bin/bash
#SBATCH --job-name=scRNAseq
#SBATCH --output=differentialabundance_%j.log
#SBATCH --error=differentialabundance_%j.err
#SBATCH --time=18:00:00  # 18 hours
#SBATCH --cpus-per-task=40  # Adjusted the number of CPUs
#SBATCH --mem=200G  # Adjusted the memory
#SBATCH --partition=scu-cpu  # The specified partition

# Run cellranger multi
/home/sor4003/store_sor4003/software_folder/cellranger-9.0.1/bin/cellranger multi \
  --id=split_CO_DiC \
  --csv=config.csv \
  --output-dir=/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/CO_DiO_DAY100_ScRNA/20250212_AV241602_2_12_2025_A_Manoj/Samples/Manoj/2_direct_analysis_cellranger_multi/outs

### config.csv

[gene-expression]
reference,/home/sor4003/store_sor4003/2a_cellranger_genome_index_nexflow/refdata-gex-GRCh38-2024-A
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
CO_DAY100,/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/CO_DiO_DAY100_ScRNA/20250212_AV241602_2_12_2025_A_Manoj/Samples/Manoj/2_direct_analysis_cellranger_multi,Gene Expression

[samples]
sample_id,ocm_barcode_ids
co_org,OB1|OB2
dic_org,OB3|OB4




########################
##########################################
ssh sor4003@scu-login02.med.cornell.edu "tar czf - --exclude='*.bam' -C /home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/CO_DiO_DAY100_ScRNA/20250212_AV241602_2_12_2025_A_Manoj/Samples/Manoj/mapping_results_with_fastq/outs/outs/per_sample_outs ." | tar xzf - -C ./


####################### march22 
#########%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^&&&&&&&&&********&&&&&&&&&&&&&
############ batch2 run 
/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/2_batch_20250313_AV241602_3_13_2025_A_Manoj/Samples/Manoj_Pool

##### need editing ## slurn submission
#!/bin/bash
#SBATCH --job-name=scRNAseq
#SBATCH --output=differentialabundance_%j.log
#SBATCH --error=differentialabundance_%j.err
#SBATCH --time=18:00:00  # 18 hours
#SBATCH --cpus-per-task=40  # Adjusted the number of CPUs
#SBATCH --mem=200G  # Adjusted the memory
#SBATCH --partition=scu-cpu  # The specified partition

# Run cellranger multi
/home/sor4003/store_sor4003/software_folder/cellranger-9.0.1/bin/cellranger multi \
  --id=split_march22 \
  --csv=config.csv \
  --output-dir=/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/2_batch_20250313_AV241602_3_13_2025_A_Manoj/Samples/Manoj_Pool


#########%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^&&&&&&&&&********&&&&&&&&&&&&& config 

[gene-expression]
reference,/home/sor4003/store_sor4003/2a_cellranger_genome_index_nexflow/refdata-gex-GRCh38-2024-A
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
CO_DAY135,/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/2_batch_20250313_AV241602_3_13_2025_A_Manoj/Samples/Manoj_Pool,Gene Expression

[samples]
sample_id,cmo_ids
co_org,OB1
dic_org,OB2
co_dio,OB3|OB4

##### excluding bam and copy all
rsync -avz --exclude='*.bam' sor4003@scu-login02.med.cornell.edu:/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/2_batch_20250313_AV241602_3_13_2025_A_Manoj/Samples/Manoj_Pool/run2_results_april8/outs/per_sample_outs .







