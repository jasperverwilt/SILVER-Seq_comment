#!/bin/bash


while getopts “d:w:n:” opt; do
  case $opt in
	d) orig_DIR=$OPTARG ;;
	w) work_DIR=$OPTARG ;;
	n) number_bases=$OPTARG ;;
  esac
done

DIR=$(pwd)

i=0


find $orig_DIR  -maxdepth 1 -name "PFP*"  -printf "%f\n" | sort -n > listsamples.txt

for sample in $(cat listsamples.txt)
do
	filename=${work_DIR}/${sample}_preprocess.sh
	echo "#!/bin/bash" > $filename
	echo "#PBS -N ${sample}_preprocess" >> $filename
	echo "#PBS -m ae" >> $filename
	echo "#PBS -l nodes=1:ppn=5,walltime=1:00:00,vmem=10gb" >> $filename
	echo "" >> $filename

	echo "cd ${orig_DIR}/${sample}" >> $filename
	echo "module purge" >> $filename
	echo "ml cutadapt/1.18-intel-2018b-Python-3.6.6" >> $filename
	echo "gunzip *.fastq.gz" >> $filename
	echo "cutadapt -u $number_bases -o ${sample}_trimmed_R1.fastq ${sample}_R1.fastq" >> $filename
	echo "cutadapt -u $number_bases -o ${sample}_trimmed_R2.fastq ${sample}_R2.fastq" >> $filename
	echo "gzip *.fastq" >> $filename
	qsub $filename
done


