#!/bin/sh
while getopts “d:o:w:” opt; do
  case $opt in
	d) orig_DIR=$OPTARG ;;
	o) out_DIR=$OPTARG ;;
	w) work_DIR=$OPTARG ;;
  esac
done


#Indexes, gtf and bed files that are used throughout the script
kal_index="/data/gent/vo/000/gvo00027/RNA_seq_pipeline/kallisto_index/Kallisto_index_hg38_91_withspikes_rDNA45S"
star_index="/data/gent/vo/000/gvo00027/RNA_seq_pipeline/STAR_index/Genome_hg38_spikes_chrIS_MTr45S_star_index"
gtf="/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes_nopatches.gtf"
exon_bed="/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/ensembl_bedregions/Homo_sapiens.GRCh38.91_exons_sorted_merged_unstranded_nopatches_tabs.bed"
intron_bed="/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/ensembl_bedregions/Homo_sapiens.GRCh38.91_introns_unstranded_nopatches_tabs.bed"
intergenic_bed="/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/ensembl_bedregions/Homo_sapiens.GRCh38.91_intergenic_unstranded_nopatches_tabs.bed"

cd $work_DIR

#Now grep searches for all normal samples, while you can change this to C to take all cancer samples
find $orig_DIR  -maxdepth 1 -name "PFP*"  -printf "%f\n" | sort -n > listsamples.txt

for sample in $(cat listsamples.txt)
do 
	echo $sample
	filename=${sample}.sh

	echo "#!/bin/sh" > $filename
	echo "#PBS -N Dedup_$sample" >> $filename
	echo "#PBS -l nodes=1:ppn=12" >> $filename
	echo "#PBS -l walltime=4:00:00" >> $filename
	echo "#PBS -l mem=45gb" >> $filename
	echo "#PBS -m abe" >> $filename
	echo "" >> $filename

	i=$sample
	DIR="${work_DIR}/$i"
	echo "mkdir $DIR" >> $filename

	echo "# Copy fastq files from original directory to scratch and unzip (if necessary)" >> $filename
	echo "cp "$orig_DIR/$i/*trimmed*" $DIR/." >> $filename
	echo "cd $DIR" >> $filename
	echo "gunzip *.fastq.gz" >> $filename
	echo "" >> $filename

	echo "# Count the number of lines per fastq to calculate the number of reads (for subsampled data), should be the same for both ends" >> $filename
	echo "wc -l *R1.fastq > total_seq_reads_subsampled_R1.txt" >> $filename
	echo "wc -l *R2.fastq > total_seq_reads_subsampled_R2.txt" >> $filename
	echo "" >> $filename

	echo "# FastQC for each sample (subsampled reads), can be combined with multiQC later on locally" >> $filename
	echo "module purge" >> $filename
	echo "module load FastQC/0.11.5-Java-1.8.0_74" >> $filename
	echo "fastqc ${i}_trimmed_R1.fastq" >> $filename
	echo "fastqc ${i}_trimmed_R2.fastq" >> $filename
	echo "mkdir FASTQC" >> $filename
	echo "mv *_fastqc* FASTQC/." >> $filename
	echo "" >> $filename

	echo "module purge" >> $filename
	echo "module load STAR/2.6.0c-intel-2018a" >> $filename
	echo "mkdir ${i}_srout" >> $filename
	echo "STAR --runThreadN 10 --outFileNamePrefix ${i}_srout/ \
	--readFilesIn ${i}_trimmed_R1.fastq ${i}_trimmed_R2.fastq  \
	--genomeDir $star_index \
	--sjdbGTFfile $gtf \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMstrandField intronMotif --outReadsUnmapped Fastx \
	--twopassMode Basic \
	--outMultimapperOrder Random \
	--outSAMmultNmax 10 \
	--outSAMprimaryFlag AllBestScore \
	--outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20" >> $filename
#Star parameters will keep the positions of aligments with best score (AllBestScore) up to 10 (outSAMmultNmax 10)

	echo "gzip *.fastq" >> $filename
	echo "echo \"STAR OK\"" >> $filename
	echo "" >> $filename

	echo "cd $DIR" >>$filename
	echo "mkdir ${i}_dedup" >> $filename
	echo "cd ${i}_dedup" >> $filename
	echo "" >> $filename

	echo "# Picard to quantify and remove PCR duplicates" >> $filename
	echo "module purge" >> $filename
	echo "ml picard/2.18.5-Java-1.8.0_162" >> $filename
	echo "java -jar \${EBROOTPICARD}/picard.jar MarkDuplicates I=$DIR/${i}_srout/Aligned.sortedByCoord.out.bam O=Aligned.sortedByCoord.picard.dedup.bam VALIDATION_STRINGENCY=SILENT M=Picard_dup.metrics REMOVE_DUPLICATES=true #remark: for MarkDuplicates the bam input file needs to be coordinate sorted" >> $filename
	echo "echo \"duplicate removal OK\"" >> $filename
	echo "" >> $filename

	echo "# Sort by name, because bamtofastq need name sorted bam" >> $filename
	echo "module purge" >> $filename
	echo "module load SAMtools/1.8-intel-2018a" >> $filename
	echo "" >> $filename

	echo "samtools sort -o Aligned.sortedByName.picard.dedup.bam -n Aligned.sortedByCoord.picard.dedup.bam" >> $filename
	echo "echo \"sorting by name OK\"" >> $filename
	echo "" >> $filename

	echo "# RSeQC: to retrieve % correct strandedness" >> $filename
	echo "module purge" >> $filename
	echo "# module load Python/2.7.14-intel-2018a" >> $filename
	echo "module load bx-python/0.8.1-intel-2018a-Python-2.7.14" >> $filename
	echo "module load RSeQC/2.6.4-intel-2018a-Python-2.7.14" >> $filename
	echo "infer_experiment.py -r /user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/ensembl_bedregions/Homo_sapiens.GRCh38.91_exons_sorted_merged2.bed -i Aligned.sortedByCoord.picard.dedup.bam > RSeQC_output_all.txt" >> $filename
	echo "" >> $filename
	echo "bam_stat.py -i Aligned.sortedByCoord.picard.dedup.bam > bam.statistics.txt" >> $filename
	echo "" >> $filename

	echo "# Convert to fastq" >> $filename
	echo "module purge" >> $filename
	echo "module load BEDTools/2.27.1-intel-2018a" >> $filename
	echo "bedtools bamtofastq -i Aligned.sortedByName.picard.dedup.bam -fq ${i}_trimmed_dedup_R1.fastq -fq2 ${i}_trimmed_dedup_R2.fastq #bam needs to be sorted by name so that pairs are together" >> $filename
	echo "# Count the number of lines per fastq to calculate the number of reads, should be the same for both ends" >> $filename
   	echo "wc -l ${i}_trimmed_dedup_R1.fastq > total_seq_reads_dedup_R1.txt" >> $filename
   	echo "wc -l ${i}_trimmed_dedup_R2.fastq > total_seq_reads_dedup_R2.txt" >> $filename
	echo "echo \"conversion of bam to fastq OK\"" >> $filename
	echo "" >> $filename

	echo "# FastQC for each sample (subsampled and duplicates removed), can be combined with multiQC later on locally" >> $filename
	echo "module purge" >> $filename
	echo "module load FastQC/0.11.5-Java-1.8.0_74" >> $filename
	echo "fastqc ${i}_trimmed_dedup_R1.fastq" >> $filename
	echo "fastqc ${i}_trimmed_dedup_R2.fastq" >> $filename
 	echo "mkdir FASTQC_dedup" >> $filename
	echo "mv *_fastqc* FASTQC_dedup/." >> $filename
	echo "" >> $filename

	echo "module purge" >> $filename
	echo "module load SAMtools/1.6-intel-2017b" >> $filename
	echo "samtools index Aligned.sortedByCoord.picard.dedup.bam" >> $filename
	echo "samtools idxstats Aligned.sortedByCoord.picard.dedup.bam > idxstats_picard.txt" >> $filename
	
	echo "samtools view -b -U reads_non_exons_picard.bam -L $exon_bed Aligned.sortedByCoord.picard.dedup.bam > reads_exons_picard.bam" >> $filename
	echo "samtools view -c reads_exons_picard.bam > statistics_picard.txt" >> $filename
	echo "samtools view -b -U reads_non_genes_picard.bam -L $intron_bed reads_non_exons_picard.bam > reads_introns_picard.bam" >> $filename
	echo "samtools view -c reads_introns_picard.bam >> statistics_picard.txt" >> $filename
	echo "samtools view -b -L  $intergenic_bed reads_non_genes_picard.bam > reads_intergenic_picard.bam" >> $filename
	echo "samtools view -c reads_intergenic_picard.bam >> statistics_picard.txt" >> $filename
	echo "mkdir regional_coverage" >> $filename
	echo "mv statistics_picard.txt reads* regional_coverage" >> $filename
	echo "cd regional_coverage" >> $filename

	echo "module purge" >> $filename
	echo "module load BEDTools/2.27.1-intel-2018a" >> $filename
	echo "intersectBed -a $gtf -b reads_introns_picard.bam  > introns_genes.txt" >> $filename
 	echo "intersectBed -a $gtf -b reads_exons_picard.bam  > exons_genes.txt" >> $filename
	echo "cd .." >> $filename

	echo "# move output folder to output directory" >> $filename
	echo "cd ${work_DIR}" >> $filename
	echo "mv ${i} $out_DIR" >> $filename
	echo "" >> $filename

	qsub $work_DIR/$filename

done 











