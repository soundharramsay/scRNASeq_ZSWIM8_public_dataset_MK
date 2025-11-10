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

##### fastq ename 
DiO_DAY184_R1.fastq.gz DiO_DAY184_S1_L001_R1_001.fastq.gz
[sor4003@scu-login02 Manoj_Pool]$ mv DiO_DAY184_R2.fastq.gz DiO_DAY184_S1_L001_R2_001.fastq.gz

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
sample_id,ocm_barcode_ids
co_org,OB1
dic_org,OB2
co_dio,OB3|OB4

##### excluding bam and copy all
rsync -avz --exclude='*.bam' sor4003@scu-login02.med.cornell.edu:/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/2_batch_20250313_AV241602_3_13_2025_A_Manoj/Samples/Manoj_Pool/run2_results_april8/outs/per_sample_outs .


####################### May 1
#########%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^&&&&&&&&&********&&&&&&&&&&&&&
############ batch3
 md5sum -c 20250424_AV241602_4_24_2025_B_Manoj.md5
DiO_DAY184_R1.fastq.gz: OK
DiO_DAY184_R2.fastq.gz: OK

 mv DiO_DAY184_R1.fastq.gz DiO_DAY184_S1_L001_R1_001.fastq.gz
 mv DiO_DAY184_R2.fastq.gz DiO_DAY184_S1_L001_R2_001.fastq.gz 


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
  --id=split_may1 \
  --csv=config.csv \
  --output-dir=/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/3_20250424_AV241602_4_24_2025_B_Manoj_run3/Samples/Manoj_Pool/

  #########%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^&&&&&&&&&********&&&&&&&&&&&&& config 

[gene-expression]
reference,/home/sor4003/store_sor4003/2a_cellranger_genome_index_nexflow/refdata-gex-GRCh38-2024-A
create-bam,true

[libraries]  # change the fastq id 
fastq_id,fastqs,feature_types
DiO_DAY184,/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/3_20250424_AV241602_4_24_2025_B_Manoj_run3/Samples/Manoj_Pool/,Gene Expression

[samples]
sample_id,ocm_barcode_ids
DIO_c184,OB1
CO_DIO_184,OB2|OB3|OB4

############################################# day 273 and 274 
[sor4003@scu-login01 CO_DAY273_fastq]$ md5sum -c 20250723_AV241602_7_23_2025_A_Manoj.md5 
CO_DAY273_R1.fastq.gz: OK
CO_DAY273_R2.fastq.gz: OK

rename 
[sor4003@scu-login01 CO_DAY273_fastq]$ mv CO_DAY273_R1.fastq.gz CO_DAY273_S1_L001_R1_001.fastq.gz
[sor4003@scu-login01 CO_DAY273_fastq]$ mv CO_DAY273_R2.fastq.gz CO_DAY273_S1_L001_R2_001.fastq.gz


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
  --id=split_day273 \
  --csv=config.csv \
  --output-dir=/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/4_20250723_AV241602_7_23_2025_A_Manoj/Samples/Manoj/CO_DAY273_fastq

[gene-expression]
reference,/home/sor4003/store_sor4003/2a_cellranger_genome_index_nexflow/custom_genome_with_gfp_mcherry/human
create-bam,true

[libraries]  # change the fastq id 
fastq_id,fastqs,feature_types
CO_DAY273,/home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/1_UCSF_MK_et_al/3_20250424_AV241602_4_24_2025_B_Manoj_run3/Samples/Manoj_Pool/,Gene Expression

[samples]
sample_id,ocm_barcode_ids
CO_DAY273,OB1
DIO_DAY273,OB2
CO-DIO_DAY273,OB3|OB4

####################### GFP/ RFP mapping 

>GFP_pENN
atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa
>mcherry
atggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcggcatggacgagctgtacaagtaa

awk '/^>/ {if (seqlen){print seqlen}; print; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' GFP_mcherry.fa
>gfp
720
>mcherry
711

echo -e 'gfp\tunknown\texon\t1\t720\t.\t+\t.\tgene_id "gfp"; transcript_id "gfp"; gene_name "gfp"; gene_biotype "protein_coding";' > GFP_mCherry.gtf

echo -e 'mCherry\tunknown\texon\t1\t711\t.\t+\t.\tgene_id "mCherry"; transcript_id "mCherry"; gene_name "mCherry"; gene_biotype "protein_coding";' >> GFP_mCherry.gtf

####
cat GFP_mcherry.fa >> genome_mcherry_gfp.fa 
####
 grep ">" genome_mcherry_gfp.fa 
 >KI270392.1 KI270392.1
>KI270394.1 KI270394.1
>gfp
>mcherry
####
 cat GFP_mCherry.gtf >> genes.gtf

 ####
 #!/bin/bash
#SBATCH --job-name=scRNAseq
#SBATCH --output=differentialabundance_%j.log
#SBATCH --error=differentialabundance_%j.err
#SBATCH --time=18:00:00  # 18 hours
#SBATCH --cpus-per-task=40  # Adjusted the number of CPUs
#SBATCH --mem=200G  # Adjusted the memory
#SBATCH --partition=scu-cpu  # The specified partition

# Run cellranger multi
/home/sor4003/store_sor4003/software_folder/cellranger-9.0.1/bin/cellranger mkref --genome=human --fasta=genome_mcherry_gfp.fa --genes=genes_mcherry_gfp.gtf  



#############
#################
################################oct4 2025 bulk 

md5sum check 
md5sum -c 20250819_AV241602_081925Manoj.md5 
CO-DIO_273_R1_R1.fastq.gz: OK------------------file transfer error rectified 
CO-DIO_273_R1_R2.fastq.gz: OK
CO-DIO_273_R2_R1.fastq.gz: OK
CO-DIO_273_R2_R2.fastq.gz: OK
CO-DIO_273_R3_R1.fastq.gz: OK
CO-DIO_273_R3_R2.fastq.gz: OK
hC0_273_R1_R1.fastq.gz: OK
hC0_273_R1_R2.fastq.gz: OK
hC0_273_R2_R1.fastq.gz: OK
hC0_273_R2_R2.fastq.gz: OK
hC0_273_R3_R1.fastq.gz: OK
hC0_273_R3_R2.fastq.gz: OK
hCO_100_R1_R1.fastq.gz: OK
hCO_100_R1_R2.fastq.gz: OK
hCO_100_R2_R1.fastq.gz: OK
hCO_100_R2_R2.fastq.gz: OK
hCO_100_R3_R1.fastq.gz: OK
hCO_100_R3_R2.fastq.gz: OK
hDIO_100_R1_R1.fastq.gz: OK
hDIO_100_R1_R2.fastq.gz: OK
hDIO_100_R2_R1.fastq.gz: OK
hDIO_100_R2_R2.fastq.gz: OK
hDIO_100_R3_R1.fastq.gz: OK
hDIO_100_R3_R2.fastq.gz: OK
hDIO_273_R1_R1.fastq.gz: OK
hDIO_273_R1_R2.fastq.gz: OK
hDIO_273_R2_R1.fastq.gz: OK
hDIO_273_R2_R2.fastq.gz: OK
hDIO_273_R3_R1.fastq.gz: OK
hDIO_273_R3_R2.fastq.gz: OK
############################# 15 samples  >>> renamed 
sample sheet 
CO-DIO_273_rep1,./CO-DIO_273_rep1_R1.fastq.gz,./CO-DIO_273_rep1_R2.fastq.gz,auto
CO-DIO_273_rep2,./CO-DIO_273_rep2_R1.fastq.gz,./CO-DIO_273_rep2_R2.fastq.gz,auto
CO-DIO_273_rep3,./CO-DIO_273_rep3_R1.fastq.gz,./CO-DIO_273_rep3_R2.fastq.gz,auto
hC0_273_rep1,./hC0_273_rep1_R1.fastq.gz,./hC0_273_rep1_R2.fastq.gz,auto
hC0_273_rep2,./hC0_273_rep2_R1.fastq.gz,./hC0_273_rep2_R2.fastq.gz,auto
hC0_273_rep3,./hC0_273_rep3_R1.fastq.gz,./hC0_273_rep3_R2.fastq.gz,auto
hCO_100_rep1,./hCO_100_rep1_R1.fastq.gz,./hCO_100_rep1_R2.fastq.gz,auto
hCO_100_rep2,./hCO_100_rep2_R1.fastq.gz,./hCO_100_rep2_R2.fastq.gz,auto
hCO_100_rep3,./hCO_100_rep3_R1.fastq.gz,./hCO_100_rep3_R2.fastq.gz,auto
hDIO_100_rep1,./hDIO_100_rep1_R1.fastq.gz,./hDIO_100_rep1_R2.fastq.gz,auto
hDIO_100_rep2,./hDIO_100_rep2_R1.fastq.gz,./hDIO_100_rep2_R2.fastq.gz,auto
hDIO_100_rep3,./hDIO_100_rep3_R1.fastq.gz,./hDIO_100_rep3_R2.fastq.gz,auto
hDIO_273_rep1,./hDIO_273_rep1_R1.fastq.gz,./hDIO_273_rep1_R2.fastq.gz,auto
hDIO_273_rep2,./hDIO_273_rep2_R1.fastq.gz,./hDIO_273_rep2_R2.fastq.gz,auto
hDIO_273_rep3,./hDIO_273_rep3_R1.fastq.gz,./hDIO_273_rep3_R2.fastq.gz,auto




##### slurm submission 
#!/bin/bash
#SBATCH --job-name=differentialabundance
#SBATCH --output=differentialabundance_%j.log
#SBATCH --error=differentialabundance_%j.err
#SBATCH --time=18:00:00          # 18 hours
#SBATCH --cpus-per-task=40       # Number of CPUs
#SBATCH --mem=200G               # Memory
#SBATCH --partition=scu-cpu      # Partition

# === Activate Conda environment ===
source /home/sor4003/anaconda3/etc/profile.d/conda.sh
conda activate env_nf

# === Nextflow memory settings (for JVM) ===
export NXF_OPTS='-Xms1g -Xmx4g'

# === Run nf-core/rnaseq pipeline ===
NXF_VER=24.04.2 nextflow run nf-core/rnaseq \
    -r 3.14.0 \
    --input /home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/2_UCSF_MK_et_al_bulk/Samples/Manoj/sample_sheet.csv \
    --igenomes_base /home/sor4003/store_sor4003/2_star_genome_index_nexflow/genome_GRch38_genecode_45_gtf_build_May2024_manually_added_cdr1as \
    -profile singularity \
    --fasta /home/sor4003/store_sor4003/2_star_genome_index_nexflow/genome_GRch38_genecode_45_gtf_build_May2024_manually_added_cdr1as/GRCh38.primary_assembly.genome.fa \
    --gtf /home/sor4003/store_sor4003/2_star_genome_index_nexflow/genome_GRch38_genecode_45_gtf_build_May2024_manually_added_cdr1as/gencode.v46.basic.annotation.gtf \
    --outdir /home/sor4003/store_sor4003/RNAseq_results_fastq/public_datasets/2_UCSF_MK_et_al_bulk/Samples/Manoj \
    -resume \
    --save_reference \
    --featurecounts_group_type "gene_type" \
    --max_memory '200.GB' \
    --max_cpus 40 \
    --max_time '30.h'


################## Differential expression 

CO-DIO_273_rep1
CO-DIO_273_rep2
CO-DIO_273_rep3

hC0_273_rep1
hC0_273_rep2
hC0_273_rep3

hCO_100_rep1
hCO_100_rep2
hCO_100_rep3

hDIO_100_rep1
hDIO_100_rep2
hDIO_100_rep3

hDIO_273_rep1
hDIO_273_rep2
hDIO_273_rep3

CO_DAY100 vs CO_DAY273 vs CO-DiO_DAY273
DIO_DAY100 vs DIO_DAY273 vs CO-DiO_DAY273











