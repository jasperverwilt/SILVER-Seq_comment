#!/bin/sh
#PBS -N kallisto_check
#PBS -l nodes=1:ppn=3
#PBS -l walltime=4:00:00
#PBS -l mem=45gb
#PBS -m abe

kal_index="/data/gent/vo/000/gvo00027/RNA_seq_pipeline/kallisto_index/Kallisto_index_hg38_91_withspikes_rDNA45S"
gtf="/data/gent/vo/000/gvo00027/RNA_seq_pipeline/ensembl_transcriptomes/Homo_sapiens.GRCh38.91_withspikes_nopatches.gtf"
work_DIR="/data/gent/vo/000/gvo00027/vsc42458/kallisto_output"
orig_DIR="/data/gent/vo/000/gvo00027/vsc42458/output_umi"

cd $work_DIR

#Now grep searches for all normal samples, while you can change this to C to take all cancer samples
find $orig_DIR  -maxdepth 1 -name "SRR*"  -printf "%f\n" | sort -n > listsamples.txt

for sample in $(cat listsamples.txt)
do 
	echo $sample
	filename=${sample}.sh

	echo "#!/bin/sh" > $filename
	echo "#PBS -N kallisto_check" >> $filename
	echo "#PBS -l nodes=1:ppn=3" >> $filename
	echo "#PBS -l walltime=4:00:00" >> $filename
	echo "#PBS -l mem=45gb" >> $filename
	echo "#PBS -m abe" >> $filename


	echo "cd /data/gent/vo/000/gvo00027/vsc42458/kallisto_output" >> $filename
	echo "mkdir ${sample}" >> $filename
	echo "cd ${sample}" >> $filename
	echo "" >> $filename

	echo "module purge" >> $filename
	echo "module load SAMtools/1.8-intel-2018a" >> $filename
	echo "samtools view -h ${orig_DIR}/${sample}/${sample}_dedup/Aligned.sortedByCoord.umi.dedup.bam | awk -v OFS="\t" '\$0 ~ /^@/{print \$0;next;} \$6 ~ /N/' > splice_reads_SRR.sam" >> $filename
	echo "samtools view -S -b splice_reads_SRR.sam > splice_reads.bam" >> $filename

	echo "samtools sort -n -o splice_reads_sorted.bam splice_reads.bam" >> $filename
	echo "samtools sort -n -o Aligned.sortedByName.umi.dedup.bam ${orig_DIR}/${sample}/${sample}_dedup/Aligned.sortedByCoord.umi.dedup.bam" >> $filename

	#echo "module purge" >> $filename
	#echo "module load HTSeq/0.11.0-foss-2018b-Python-2.7.15" >> $filename 
	#echo "htseq-count --format bam --order name \
	#--nonunique none \
	#--stranded no \
	#Aligned.sortedByName.picard.dedup.bam \
	#$gtf > htseq_counts_total_SRR.txt" >> $filename

	#echo "htseq-count --format bam --order name \
	#--nonunique none \
	#--stranded no \
	#splice_reads_sorted.bam \
	#$gtf > htseq_counts_total_SRR_splice.txt" >> $filename

	echo "module purge" >> $filename
	echo "ml BEDTools/2.27.1-intel-2018a" >> $filename
	echo "bedtools bamtofastq -i splice_reads_sorted.bam -fq SILVER_splice.fastq" >> $filename
	echo "bedtools bamtofastq -i Aligned.sortedByName.picard.dedup.bam -fq SILVER.fastq" >> $filename

	echo "module purge" >> $filename
	echo "ml kallisto/0.44.0-intel-2018a" >> $filename
	echo "kallisto quant -t 10 -l 260 -s 30 -o klout_splice SILVER_splice.fastq" >> $filename
	echo "kallisto quant -t 10 -l 260 -s 30 -i $kal_index -o klout SILVER.fastq" >> $filename

	echo "module purge" >> $filename
	echo "ml R/3.6.0-intel-2019a" >> $filename
	echo "Rscript ../counter_kal.R >> ../kallisto.txt" >> $filename
	echo "echo ","  >> ../kallisto.txt" >> $filename
	qsub $filename
done
