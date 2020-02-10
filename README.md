# SILVER-Seq letter
We analysed publicly available data from ["Extracellular RNA in a single droplet of human serum reflects physiologic and disease states"](https://www.pnas.org/content/116/38/19200.short) for DNA contamination by calculating several metrics.

The following metrics have been evaluated: 
* Strandedness
* Coverage by genomic region
* Reads mapping to spliced sequences
* Inference of DNA copy number profile

## Methods
Raw data was downloaded from the NCBI SRA database (SRA Study SRP198979) and all the preprocessing steps were performed in accordance to the sample report on the Genemo's [SILVER-Seq](https://genemo.com/services/silver-seq/) website, unless reported differently. The UMI reads of both serum and normal samples were merged with the corresponding raw reads. The package ```umitools 1.0.0``` [[1]](#1) was then used to extract the UMIs from the reads. Next, the first five bases of each read were trimmed using ```Cutadapt 1.18``` [[2]](#2). The reads were mapped to a h38 reference genome using ```STAR 2.6.0c``` and ```umitools 1.0.0``` [[1]](#1) was used to remove PCR-duplicates. The infer_expiment.py function from ```RSeQC 2.6.4``` was used to determine strandedness of the reads. ```Samtools 1.8``` was used to determine the spliced and non-spliced nature of the reads and to calculate the fraction of reads mapping to exonic, intronic and intergenic regions. After processing the raw FASTQ files with ```Cutadapt 1.18```, ```bwa-mem 0.7.17``` and ```Picard 2.1.1```, copy number profiles were generated using ```WisecondorX 1.1.15```. The results obtained after re-analyzing SILVER-Seq data were compared with our published results on total RNA-seq data from platelet-free plasma. 

### References
<a id="1">[1]</a> 
Smith, T., Heger, A., & Sudbery, I. (2017). 
UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. 
Genome research, 27(3), 491-499.

<a id="2">[2]</a> 
Martin, M. (2011).
Cutadapt removes adapter sequences from high-throughput sequencing reads.
EMBnet. journal, 17(1), 10-12.
