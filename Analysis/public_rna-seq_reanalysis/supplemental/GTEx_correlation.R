# GTEx correlation analysis

library(tidyverse)
library(readxl)
library(gprofiler2)
library(data.table)
GTEx_table <- data.table::fread('~/Documents/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct') %>%
  separate(Name, into = c('ENSG', 'ver'))
mainData
conv_GS_targs <- gconvert(
  query = mainData$Gene,
  organism = "hsapiens",
  target = "ENSG",
  numeric_ns = "",
  mthreshold = 1,
  filter_na = FALSE
) # SOX10 doesn't exist in GTEx table; exclude from analysis
conv_GS_paras <- gconvert(
  query = mainData$Paralog,
  organism = "hsapiens",
  target = "ENSG",
  numeric_ns = "",
  mthreshold = 1,
  filter_na = FALSE
)

mainData_GTEx <- mainData
mainData_GTEx$ENSG_targ <- conv_GS_targs$target
mainData_GTEx$ENSG_para <- conv_GS_paras$target
pairdat <- tibble(ENSG_targ = character(),
                  Description_targ = character(),
                  tissue = character(),
                  logTPMp1_targ = numeric(),
                  ENSG_para = character(),
                  Description_para = character(),
                  logTPMp1_para = numeric())
paircor_sum <- tibble(Gene = character(),
                          Paralog = character(),
                          FC2 = numeric(),
                          padj = numeric(),
                          ENSG_targ = character(),
                          ENSG_para = character(),
                          cor_spearman = numeric())
for (pair in 1:nrow(mainData_GTEx)) {
  
  tmp_maindat <- mainData_GTEx[pair,]
  if(tmp_maindat$Gene == 'SOX10') {next}
  tmp_GTEx_targ <- GTEx_table %>%
    filter(ENSG == tmp_maindat$ENSG_targ) %>%
    pivot_longer("Adipose - Subcutaneous":"Whole Blood", names_to = 'tissue', values_to = 'TPM') %>%
    mutate(logTPMp1_targ = log(TPM+1)) %>%
    dplyr::select(-c(ver, TPM)) %>%
    dplyr::rename(ENSG_targ = ENSG,
                  Description_targ = Description)
  tmp_GTEx_para <- GTEx_table %>%
    filter(ENSG == tmp_maindat$ENSG_para) %>%
    pivot_longer("Adipose - Subcutaneous":"Whole Blood", names_to = 'tissue', values_to = 'TPM') %>%
    mutate(logTPMp1_para = log(TPM+1)) %>%
    dplyr::select(-c(ver, TPM)) %>%
    dplyr::rename(ENSG_para = ENSG,
                  Description_para = Description)
  
  pairdat_tmp <- inner_join(tmp_GTEx_targ, tmp_GTEx_para, by = 'tissue')
  
  pairdat %<>% bind_rows(pairdat_tmp)
  
  cor_tmp <- cor(pairdat_tmp$logTPMp1_targ, pairdat_tmp$logTPMp1_para, method = 'spearman')
  
  paircor_sum_tmp <- tmp_maindat %>%
    dplyr::select(Gene, Paralog, FC2, padj) %>%
    mutate(ENSG_targ = unique(pairdat_tmp$ENSG_targ),
           ENSG_para = unique(pairdat_tmp$ENSG_para),
           cor_spearman = cor_tmp)
  
  paircor_sum %<>% bind_rows(paircor_sum_tmp)
  
}
ggplot(paircor_sum, aes(cor_spearman, FC2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
cor.test(paircor_sum$FC2, paircor_sum$cor_spearman, method = 'spearman')

write.csv(paircor_sum, file = '~/code/GitHub/grn_nitc/resub1/Analysis/public_rna-seq_reanalysis/supplemental/GTEx_pariwisecor_FC2.csv', quote = F, row.names = F)
