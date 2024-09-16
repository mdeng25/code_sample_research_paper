options(stringsAsFactors = FALSE)
theme_set(figure_theme)

require(tidyverse)

setwd("/Users/Mark/Desktop/Research/ScreenPaper/Ribosome profiling/Jan2020_Expt1/")

source("Lance code/CodonInfo.R")

## Import data

codon.counts.1 <- tbl_df(read.csv(file="Lance code/codon_counts_20200401/NI810_codon_counts.tsv", sep="\t")) %>% mutate(sample = "B_D")
codon.counts.2 <- tbl_df(read.csv(file="Lance code/codon_counts_20200401/NI811_codon_counts.tsv", sep="\t")) %>% mutate(sample = "B_L1")
codon.counts.3 <- tbl_df(read.csv(file="Lance code/codon_counts_20200401/NI812_codon_counts.tsv", sep="\t")) %>% mutate(sample = "B_L24")
codon.counts.4 <- tbl_df(read.csv(file="Lance code/codon_counts_20200401/NI813_codon_counts.tsv", sep="\t")) %>% mutate(sample = "KO_D")
codon.counts.5 <- tbl_df(read.csv(file="Lance code/codon_counts_20200401/NI814_codon_counts.tsv", sep="\t")) %>% mutate(sample = "KO_L1")
codon.counts.6 <- tbl_df(read.csv(file="Lance code/codon_counts_20200401/NI815_codon_counts.tsv", sep="\t")) %>% mutate(sample = "KO_L24")

## Combine data

codon.counts.epa.all <- rbind(codon.counts.1, codon.counts.2, codon.counts.3, codon.counts.4, codon.counts.5, codon.counts.6)

codon.counts.epa.geneids <- tbl_df(data.frame(gene_id = unique(codon.counts.epa.all$gene_id)))

## Convert Ensembl IDs to Gene Symbols

require("org.Mm.eg.db")
key_field = "ENSEMBL"
cols <- c("ALIAS","SYMBOL", "ENTREZID", "GENENAME")
last <- function(x){x[[length(x)]]}
last.fancy <- function (x) {
  ret <- x[[length(x)]]
  if(length(x) > 1 & substr(ret,1,3) == "LOC") {ret <- x[[length(x)-1]]}
  if(length(x) > 2 & substr(ret,1,3) == "LOC") {ret <- x[[length(x)-2]]}
  ret
}
for (col in cols) {
  codon.counts.epa.geneids[,col] <- mapIds(org.Mm.eg.db, keys=as.character(codon.counts.epa.geneids$gene_id), 
                                           column=col, keytype=key_field, multiVals=last.fancy)
}

codon.counts.epa.all.mg <- tbl_df(merge(codon.counts.epa.geneids, codon.counts.epa.all, by="gene_id"))

## AFTER RUNNING THIS CODE ONCE, RUN THE FOLLOWING LINE. NEXT TIME, START FROM THE LINE BELOW AND AVOID RUNNING THE ABOVE CODE
# write.csv(codon.counts.epa.all.mg, file="codon.counts.epa.all.mg.csv") # Note this file will be huge ~2.5 GB
# same.data.imported <- tbl_df(read.csv(file="codon.counts.epa.all.mg.csv")) %>% dplyr::select(-X)
# codon.counts.epa.all.mg <- same.data.imported

# codon.counts.epa.all.mg %>% filter(gene_id == "ENSMUSG00000028495")

## Rename columns for convenience

codon.counts.epa.rename <- codon.counts.epa.all.mg %>% mutate(asite_p1 = esite_position_1_count, asite_p2 = esite_position_2_count, asite_p3 = esite_position_3_count,
                                                           psite_p1 = psite_position_1_count, psite_p2 = psite_position_2_count, psite_p3 = psite_position_3_count,
                                                           esite_p1 = asite_position_1_count, esite_p2 = asite_position_2_count, esite_p3 = asite_position_3_count) %>% 
  mutate(codon_seq = substr(codon_seq, 7, 9)) %>% 
  dplyr::select(gene_id, ENTREZID, ALIAS, GENENAME, gene_length_bp, gene_length_codons, sample, codon_seq, codon_index, codon_count_sum, asite_p1:esite_p3)
codon.counts.epa.rename

codon.counts.epa.rename.upper <- codon.counts.epa.rename %>% filter(codon_seq == toupper(codon_seq))

codon.counts.basic.stats <- codon.counts.epa.rename.upper %>% group_by(gene_id, sample) %>% 
  mutate(sum = sum(codon_count_sum)/3, p1_sum = sum(asite_p1), p2_sum = sum(asite_p2), p3_sum = sum(asite_p3))
codon.counts.filter <- codon.counts.basic.stats %>% filter(!is.na(ALIAS), substr(ALIAS,1,5) != "Snrnp") 
codon.counts.filter.1000 <- codon.counts.filter %>% group_by(gene_id) %>% filter(sum(sum) > 1000)

codon.counts.filter %>% filter(ALIAS == "Actb")
# as.data.frame(codon.counts.filter %>% filter(ALIAS == "Capzb", sample == "KO_L1") %>% dplyr::select(-GENENAME))[1:10,]
# as.data.frame(codon.counts.filter %>% filter(ALIAS == "Capzb", sample == "KO_L1") %>% filter(codon_seq %in% c("CTC","CTT")) %>% dplyr::select(-GENENAME))

codon.counts.filter.Mult3 <- codon.counts.filter %>% group_by(gene_id, gene_length_bp) %>% filter(gene_length_bp%%3 == 0)
codon.counts.filter.NotMult3 <- codon.counts.filter %>% group_by(gene_id, gene_length_bp) %>% filter(gene_length_bp%%3 != 0)
Mult3.genes <- unique(codon.counts.filter.Mult3$gene_id)

## Import codon frequency data
## BELOW IS CODE THAT FIGURES OUT WHAT THE CODON FREQUENCIES ARE IN THE GENES THAT PASSED QUALITY CONTROL
## not how many ribosomes are bound, just how many codons of each kind appear in these genes...

# Use codon tallies to make a table of codon frequencies
codon.tallies.correct <- tbl_df(read.csv(file="Lance code/ByGene.Codon_frequencies", sep="\t"))
codon.tallies.correct.Mult3 <- codon.tallies.correct %>% filter(gene_id %in% Mult3.genes)
codon.tallies.correct.gat <- codon.tallies.correct.Mult3 %>% gather(codon_seq, count, -gene_id, -gene_length_bp, -gene_length_codons)
codon.tallies.correct.mut <- codon.tallies.correct.gat %>% group_by(gene_id) %>% 
  mutate(gene.count = sum(count), codon.freq = count/gene.count) 

codon.tallies.allgenes <- codon.tallies.correct.gat %>% group_by(codon_seq) %>% 
  summarize(codon.sum = sum(count)) %>% ungroup %>% mutate(codon.freq = codon.sum/sum(codon.sum))
# as.data.frame(codon.tallies.allgenes)
codon.tallies.allgenes.sel <- codon.tallies.allgenes %>% dplyr::select(-codon.sum)

