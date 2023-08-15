## Run once to generate Resources/ contents
# See README for file sources

library(GenomicFeatures)

# Get "gene lengths", (i.e., union of all Refseq transcript exons/UTRs per gene) from hg38 genes
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ 
# file: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz on 2022/05/05
txdb <- makeTxDbFromGFF("../hg38.ncbiRefSeq.gtf", format="gtf")
lengthsPergeneid <- sum(width(IRanges::reduce(exonsBy(txdb, by = "gene")))) # a named int(), in bases
lengthtbl<-as_tibble(list(
  gene_id = names(lengthsPergeneid),
  length = lengthsPergeneid
))
write.csv(lengthtbl, file = '~/code/grn_nitc/Resources/hg38.ncbiRefSeq.txLengthPerGene.csv', row.names = F)
