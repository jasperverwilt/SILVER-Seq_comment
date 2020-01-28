#!/bin/bash


while getopts “d:s:w:p:n:” opt; do
  case $opt in
	d) orig_DIR=$OPTARG ;;
	s) SRR_DIR=$OPTARG ;;
	w) work_DIR=$OPTARG ;;
	p) pattern=$OPTARG ;; 
	n) number_bases=$OPTARG ;;
  esac
done

DIR=$(pwd)

i=0

for sample in $(cat ${DIR}/SRR_Acc_List.txt)
do
	filename=${work_DIR}/${sample}_download.sh
	i=$((i+1))
	umi=$(sed "${i}q;d" ${DIR}/SRR_Acc_List_UMI.txt)
	echo "#!/bin/bash" > $filename
	echo "#PBS -N download_${umi}" >> $filename
	echo "#PBS -m ae" >> $filename
	echo "#PBS -l nodes=1:ppn=3,walltime=2:00:00,vmem=10gb" >> $filename
	echo "" >> $filename

	echo "cd $orig_DIR" >> $filename
	echo "module purge" >> $filename
	echo "ml SRA-Toolkit/2.9.6-1-centos_linux64" >> $filename
	echo "" >> $filename
	
	echo "mkdir ${sample}" >> $filename
	echo "prefetch -v $sample" >> $filename
	echo "prefetch -v $umi" >> $filename
	echo "cd ${sample}" >> $filename
	echo "fastq-dump ${SRR_DIR}/${umi}.sra" >> $filename
	echo "fastq-dump ${SRR_DIR}/${sample}.sra" >> $filename
	echo "mv ${umi}.fastq ${sample}_umi.fastq" >> $filename

	echo "ml Python/3.6.6-intel-2018b" >> $filename
	echo "python3 ${DIR}/move_UMI.py ${orig_DIR}/${sample}/${sample}.fastq ${orig_DIR}/${sample}/${sample}_umi.fastq > ${sample}_new.fastq" >> $filename
	echo "gzip ${sample}.fastq" >> $filename
	echo "gzip ${sample}_umi.fastq" >> $filename

	echo "module purge" >> $filename
	echo "ml UMI-tools/1.0.0-intel-2018b-Python-3.6.6" >> $filename
	echo "" >> $filename

	echo "gunzip *_new.fastq.gz" >> $filename
	echo "umi_tools extract --stdin=${sample}_new.fastq --bc-pattern=$pattern --log=processed.log --stdout ${sample}_prcs.fastq.gz" >> $filename
	echo "gzip ${sample}_new.fastq" >> $filename

	echo "module purge" >> $filename
	echo "ml cutadapt/1.18-intel-2018b-Python-3.6.6" >> $filename
	echo "" >> $filename

	echo "gunzip *_prcs.fastq.gz" >> $filename
	echo "cutadapt -u $number_bases -o ${sample}_trimmed.fastq ${sample}_prcs.fastq" >> $filename
	echo "gzip *.fastq" >> $filename
	qsub $filename
done


