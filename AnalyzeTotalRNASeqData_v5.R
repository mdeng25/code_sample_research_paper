setwd("/Users/Mark/Desktop/Research/ScreenPaper/Ribosome profiling/Jan2020_Expt1/")

require(tidyverse)

## IMPORT DATA AND MERGE

dat1 <- tbl_df(read.table(file="total rnaseq data/counts_1.txt", header=TRUE))
dat2 <- tbl_df(read.table(file="total rnaseq data/counts_2.txt", header=TRUE))
dat3 <- tbl_df(read.table(file="total rnaseq data/counts_3.txt", header=TRUE))
dat4 <- tbl_df(read.table(file="total rnaseq data/counts_4.txt", header=TRUE))
dat5 <- tbl_df(read.table(file="total rnaseq data/counts_5.txt", header=TRUE))
dat6 <- tbl_df(read.table(file="total rnaseq data/counts_6.txt", header=TRUE))

dat1c <- dat1 %>% mutate(WT_D_m = AlignReads2.Aligned_1.Aligned.sortedByCoord.out.bam) %>% 
  dplyr::select(-AlignReads2.Aligned_1.Aligned.sortedByCoord.out.bam)
dat2c <- dat2 %>% mutate(WT_L1_m = AlignReads2.Aligned_2.Aligned.sortedByCoord.out.bam) %>% 
  dplyr::select(-AlignReads2.Aligned_2.Aligned.sortedByCoord.out.bam)
dat3c <- dat3 %>% mutate(WT_L24_m = AlignReads2.Aligned_3.Aligned.sortedByCoord.out.bam) %>% 
  dplyr::select(-AlignReads2.Aligned_3.Aligned.sortedByCoord.out.bam)
dat4c <- dat4 %>% mutate(KO_D_m = AlignReads2.Aligned_4.Aligned.sortedByCoord.out.bam) %>% 
  dplyr::select(-AlignReads2.Aligned_4.Aligned.sortedByCoord.out.bam)
dat5c <- dat5 %>% mutate(KO_L1_m = AlignReads2.Aligned_5.Aligned.sortedByCoord.out.bam) %>% 
  dplyr::select(-AlignReads2.Aligned_5.Aligned.sortedByCoord.out.bam)
dat6c <- dat6 %>% mutate(KO_L24_m = AlignReads2.Aligned_6.Aligned.sortedByCoord.out.bam) %>% 
  dplyr::select(-AlignReads2.Aligned_6.Aligned.sortedByCoord.out.bam)

dat.merged.12 <- tbl_df(merge(dat1c, dat2c, by=c("Geneid","Chr","Start","End","Strand","Length")))
dat.merged.34 <- tbl_df(merge(dat3c, dat4c, by=c("Geneid","Chr","Start","End","Strand","Length")))
dat.merged.56 <- tbl_df(merge(dat5c, dat6c, by=c("Geneid","Chr","Start","End","Strand","Length")))

dat.merged.0123 <- tbl_df(merge(dat.merged.12, dat.merged.34, by=c("Geneid","Chr","Start","End","Strand","Length")))
dat.merged <- tbl_df(merge(dat.merged.0123, dat.merged.56, by=c("Geneid","Chr","Start","End","Strand","Length")))

# dat.merged %>% dplyr::select(-Start, -End, -Strand, -Length, -Chr) %>% arrange(-c1)

## Geneid's are ENSEMBL ID's (look up ENSEMBL), which are unique identifiers but impossible to understand
## The following chunk of code merges in classical gene names that are easier to deal with

## MERGE IN GENE NAMES, ETC

require("org.Mm.eg.db")
key_field = "ENSEMBL"
cols <- c("ALIAS", "ENTREZID", "GENENAME")
last <- function(x){x[[length(x)]]}
last.fancy <- function (x) {
  ret <- x[[length(x)]]
  if(length(x) > 1 & substr(ret,1,3) == "LOC") {ret <- x[[length(x)-1]]}
  if(length(x) > 2 & substr(ret,1,3) == "LOC") {ret <- x[[length(x)-2]]}
  ret
}
for (col in cols) {
  dat.merged[,col] <- mapIds(org.Mm.eg.db, keys=as.character(dat.merged$Geneid),
                             column=col, keytype=key_field, multiVals=last.fancy)
}

dat.colnames.mRNA <- dat.merged %>% dplyr::select(Geneid, ENTREZID, ALIAS, GENENAME, WT_D_m:KO_L24_m)

## Un-comment the following lines and running them to get a sense of the data
# dat.colnames.mRNA %>% filter(ALIAS == "Actb")
# dat.colnames.mRNA %>% arrange(-WT_D_m)

## Saving this in a csv, so that we don't have to re-run the code above if we don't want to
# write.csv(dat.colnames.mRNA, "Total_mRNA_data_20200408.csv", row.names = F)

## Saving an even simpler file, without ENSEMBL and ENTREZ ID's. Use the one above, if anything (throwing away the ID's is pointless, and they may come in handy)

dat.colnames.mRNA.mod <- dat.colnames.mRNA %>% mutate(Symbol = ALIAS, Annot = GENENAME) %>% dplyr::select(-Geneid, -ENTREZID, -ALIAS, -GENENAME)
dat.colnames.mRNA.mod <- dat.colnames.mRNA.mod %>% dplyr::select(Symbol, Annot, WT_D_m:KO_L24_m)
# write.csv(dat.colnames.mRNA.mod, "Total_mRNA_GCN2_data_20200408_simpler.csv", row.names = F)
