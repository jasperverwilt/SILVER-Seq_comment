# Pipeline for SILVER-Seq data

## Step 1: Download and preprocess
This bash scripts writes a SBATCH jobfile for each sample to your a temporary directory and submits it to the cluster. The samples you wish to download should be in a file called ```SRR_Acc_List.txt``` and the corresponding UMIs in a file called ```SRR_Acc_List_UMI.txt```. The script takes the following inputs:
* ```-d```: the directory you wish to save the downloaded and processed FASTQ files in
* ```-s```: the directory where your .sra files are stored
* ```-w```: the directory which will contain your job scripts
* ```-p```: the pattern that should be used to extract the UMIs
* ```-n```: the number of bp that should be trimmed from the reads (without UMIs)

Take a look at the following example use:
```bash
bash SILVER_preprocessing.sh -d /data/gent/vo/000/gvo00027/vsc42458/origin_test -s /data/gent/vo/000/gvo00027/vsc42458/SRR_files/sra/ -w /scratch/gent/vo/000/gvo00027/tmp_output_FeO/ -p NNNNNNNN -n 5
```

The individual jobscripts will perform the following steps: 
1. Download the raw reads and UMI FASTQ files (```SRR9094524.fastq.gz```and ```SRR9094524_umi.fastq.gz```)
2. Merges each read with its corresponding UMI (```SRR9094524_new.fastq.gz```)
3. Extracts the UMI and adds it to the header (```SRR9094524_prcs.fastq.gz```)
4. Cut the first number of bases (```SRR9094524_trimmed.fastq.gz```)

As an example the output of the script is printed for SRR9094524:
```bash
SRR9094524
├── processed.log
├── SRR9094524.fastq.gz
├── SRR9094524_new.fastq.gz
├── SRR9094524_prcs.fastq.gz
├── SRR9094524_trimmed.fastq.gz
└── SRR9094524_umi.fastq.gz
```

### UMItools and Cutadapt
UMIs were provided in separate FASTQ files, so had to be merged with the raw reads for deduplification later on. For the cancer samples the reads in the UMI FASTQ had a lenght of 16, with 1-8 being a barcode and 9-16 being the actual UMI. Therefore, the pattern that needs to be provided is "XXXXXXXXNNNNNNNN". That way only the first 9-16 bases will be recognized as UMI. This does mean however that you would have to add eight to the number of bases you want to trim from the raw reads. For the normal samples, the UMI spanned the first eight bases, so the pattern is "NNNNNNNN". We trimmed the first five bases from the actual reads, in correspondance with the workflow the authors published in a sample report on the [website of SILVER-Seq](https://genemo.com/services/silver-seq/). 

## Step 2: Generate RNA-seq QC parameters
Now that all the preprocessing is done, this script creates and runs a jobfile for each sample in the list. You can submit the this script and change the origin, working and and output directory in the code of the script itself. The script takes the following inputs:
* ```-d```: the directory containing the downloaded and preprocessed data
* ```-o```: the output directory
* ```-w```: the directory which will contain your job scripts

Take a look at the following example use:
```bash
bash SILVER_data_check.sh -d /data/gent/vo/000/gvo00027/vsc42458/origin_test -o /data/gent/vo/000/gvo00027/vsc42458/output_test -w /scratch/gent/vo/000/gvo00027/tmp_output_FeO/ 
```

The steps can be generalized as follows:
1. ```FastQC``` generation of duplicated data
2. ```STAR``` mapping
3. Deduplification with ```UMItools```
4. Calculate strandedness and splice reads with ```RSeQC```
5. ```FastQC``` generation of deduplicated data
6. Find which reads map to exonic, intronic and intergenic regions
7. Find which genes the exonic and intronic reads map to

The output for SRR9094428 will look like this:

```bash
SRR9094524
├── FASTQC
│   ├── SRR9094524_trimmed_fastqc.html
│   └── SRR9094524_trimmed_fastqc.zip
├── SRR9094524_dedup
│   ├── Aligned.sortedByCoord.umi.dedup.bam
│   ├── Aligned.sortedByCoord.umi.dedup.bam.bai
│   ├── Aligned.sortedByName.umi.dedup.bam
│   ├── bam.statistics.txt
│   ├── deduplicated_edit_distance.tsv
│   ├── deduplicated_per_umi_per_position.tsv
│   ├── deduplicated_per_umi.tsv
│   ├── FASTQC_dedup
│   │   ├── SRR9094524_dedup_fastqc.html
│   │   └── SRR9094524_dedup_fastqc.zip
│   ├── idxstats_umi.txt
│   ├── regional_coverage
│   │   ├── exons_genes.txt
│   │   ├── introns_genes.txt
│   │   ├── reads_exons_umi.bam
│   │   ├── reads_intergenic_umi.bam
│   │   ├── reads_introns_umi.bam
│   │   ├── reads_non_exons_umi.bam
│   │   ├── reads_non_genes_umi.bam
│   │   └── statistics_umi.txt
│   ├── RSeQC_output_all.txt
│   ├── SRR9094524_dedup.fastq
│   └── total_seq_reads_dedup.txt
├── SRR9094524_srout
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.sorted.bam
│   ├── Aligned.sortedByCoord.out.sorted.bam.bai
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
│   └── Unmapped.out.mate1
├── SRR9094524_trimmed.fastq.gz
└── total_seq_reads_subsampled.txt
```
