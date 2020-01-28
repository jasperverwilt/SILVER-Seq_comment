# Copy Number Profiles
1. We downloaded the fastq files from SILVER-seq and performed trimming
2. WisecondorX (downloaded from github on 14/01/2020) was run on all the breast cancer cells, after creating a panel of normal with the "normal" samples, except for SRR9094547, with the command `WisecondorX newref <files in normalFiles.txt> SILVER.GRCh38.200kbp.noNIPT.npz --cpus 3 --binsize 200000`
    * The full pipeline of WisecondorX is found in 20200114_sWGS_pipeline_SE.snakefile
3. We used a modified plotter.R to plot SRR9094547 (adjustedPlotter.R), so that the gonosomes are found on 1n, relative to the baseline (and not on the baseline)