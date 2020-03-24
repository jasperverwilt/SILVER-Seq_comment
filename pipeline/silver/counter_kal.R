library(tidyverse)
genes <- read.table("klout/abundance.tsv", header = TRUE)
genes_splice <- read.table("klout_splice/abundance.tsv", header = TRUE)

spliced_genes <- genes_splice %>%
  filter(tpm > 0) %>%
  .$target_id

cat(nrow(genes %>%
       filter(target_id %in% spliced_genes, tpm > 5)))
cat(getwd())
cat(nrow(filter(genes, tpm > 5)))
