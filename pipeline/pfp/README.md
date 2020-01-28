# Pipeline for PFP data

## Step 1: Preprocess
This bash scripts writes a SBATCH jobfile for each sample to your a temporary directory and submits it to the cluster. The script takes the following inputs:
* ```-d```: the directory you wish to save the downloaded and processed FASTQ files in
* ```-w```: the directory which will contain your job scripts
* ```-n```: the number of bp that should be trimmed from the reads (without UMIs)

Take a look at the following example use:
```bash
bash PFP_preprocessing.sh -d /data/gent/vo/000/gvo00027/vsc42458/SRR_samples/data_CEV -w /scratch/gent/vo/000/gvo00027/tmp_output_FeO/ -n 5
```
All this script does is use ```Cutadapt```to remove the first number of bases and save the new FASTQ files ending in ```_trimmed_R1.fastq.gz```and ```_trimmed_R2.fastq.gz```.

As an example the output of the script is printed for SRR9094524:
```bash
PFP_rep1_donor1
├── PFP_rep1_donor1_R1.fastq.gz
├── PFP_rep1_donor1_R2.fastq.gz
├── PFP_rep1_donor1_trimmed_R1.fastq.gz
└── PFP_rep1_donor1_trimmed_R2.fastq.gz
```

## Step 2: Generate RNA-seq QC parameters
Now that all the preprocessing is done, this script creates and runs a jobfile for each sample in the list. You can submit the this script and change the origin, working and and output directory in the code of the script itself. The script takes the following inputs:
* ```-d```: the directory containing the downloaded and preprocessed data
* ```-o```: the output directory
* ```-w```: the directory which will contain your job scripts

Take a look at the following example use:
```bash
bash PFP_data_check.sh -d /data/gent/vo/000/gvo00027/vsc42458/SRR_samples/data_CEV -o /data/gent/vo/000/gvo00027/vsc42458/output_test -w /scratch/gent/vo/000/gvo00027/tmp_output_FeO/ 
```

The steps can be generalized as follows:
1. ```FastQC``` generation of duplicated data
2. ```STAR``` mapping
3. Deduplification with ```Picard```
4. Calculate strandedness and splice reads with ```RSeQC```
5. ```FastQC``` generation of deduplicated data
6. Find which reads map to exonic, intronic and intergenic regions
7. Find which genes the exonic and intronic reads map to

The output for PFP_rep1_donor1 will look like this:

```bash
PFP_rep1_donor1
├── FASTQC
│   ├── PFP_rep1_donor1_trimmed_R1_fastqc.html
│   ├── PFP_rep1_donor1_trimmed_R1_fastqc.zip
│   ├── PFP_rep1_donor1_trimmed_R2_fastqc.html
│   └── PFP_rep1_donor1_trimmed_R2_fastqc.zip
├── PFP_rep1_donor1_dedup
│   ├── Aligned.sortedByCoord.picard.dedup.bam
│   ├── Aligned.sortedByCoord.picard.dedup.bam.bai
│   ├── Aligned.sortedByName.picard.dedup.bam
│   ├── bam.statistics.txt
│   ├── FASTQC_dedup
│   │   ├── PFP_rep1_donor1_trimmed_dedup_R1_fastqc.html
│   │   ├── PFP_rep1_donor1_trimmed_dedup_R1_fastqc.zip
│   │   ├── PFP_rep1_donor1_trimmed_dedup_R2_fastqc.html
│   │   └── PFP_rep1_donor1_trimmed_dedup_R2_fastqc.zip
│   ├── idxstats_picard.txt
│   ├── PFP_rep1_donor1_trimmed_dedup_R1.fastq
│   ├── PFP_rep1_donor1_trimmed_dedup_R2.fastq
│   ├── Picard_dup.metrics
│   ├── regional_coverage
│   │   ├── exons_genes.txt
│   │   ├── introns_genes.txt
│   │   ├── reads_exons_picard.bam
│   │   ├── reads_intergenic_picard.bam
│   │   ├── reads_introns_picard.bam
│   │   ├── reads_non_exons_picard.bam
│   │   ├── reads_non_genes_picard.bam
│   │   └── statistics_picard.txt
│   ├── RSeQC_output_all.txt
│   ├── total_seq_reads_dedup_R1.txt
│   └── total_seq_reads_dedup_R2.txt
├── PFP_rep1_donor1_srout
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   ├── _STARgenome
│   │   ├── exonGeTrInfo.tab
│   │   ├── exonInfo.tab
│   │   ├── geneInfo.tab
│   │   ├── sjdbInfo.txt
│   │   ├── sjdbList.fromGTF.out.tab
│   │   ├── sjdbList.out.tab
│   │   └── transcriptInfo.tab
│   ├── _STARpass1
│   │   ├── Log.final.out
│   │   └── SJ.out.tab
│   ├── Unmapped.out.mate1
│   └── Unmapped.out.mate2
├── PFP_rep1_donor1_trimmed_R1.fastq.gz
├── PFP_rep1_donor1_trimmed_R2.fastq.gz
├── total_seq_reads_subsampled_R1.txt
└── total_seq_reads_subsampled_R2.txt
```
