library(tidyverse)
library(biomaRt)
library(magrittr)
library(cvequality)

source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')

setwd('~/code/grn_nitc/rnaseq/')

datasets <- as_tibble(read.csv('annotations/dataset_metadata.csv', stringsAsFactors = F, header = T))

foldchanges_hitsonly_p_fc <- as_tibble(read.csv('de_analysis/summary_plots_drafts/foldchanges_hitsonly_p_fc.csv', stringsAsFactors = F, header = T))
foldchanges_nonhitsonly_anypara <- as_tibble(read.csv('de_analysis/summary_plots_drafts/foldchanges_nonhitsonly_anypara.csv', stringsAsFactors = F, header = T))

tflist_f_Hs <- '~/code/grn_nitc/Resources/Homo_sapiens_TF.txt'
tflist_f_Mm <- '~/code/grn_nitc/Resources/Mus_musculus_TF.txt'

tflist <- as_tibble(read.table(tflist_f_Hs, header = T, sep = '\t')) %>%
  dplyr::rename(CRISPR_target_GS = Symbol) %>%
  bind_rows(as_tibble(read.table(tflist_f_Mm, header = T, sep = '\t')) %>%
              dplyr::rename(CRISPR_target_GS = Symbol))

hits_TF_para_pairs_p_fc <- foldchanges_hitsonly_p_fc %>% filter(KO_Gene %in% tflist$CRISPR_target_GS, Paralog %in% tflist$CRISPR_target_GS)
nonhits_TF_para_pairs_anypara <- foldchanges_nonhitsonly_anypara %>% filter(KO_Gene %in% tflist$CRISPR_target_GS, Paralog %in% tflist$CRISPR_target_GS)

write.csv(hits_TF_para_pairs_p_fc, file = 'de_analysis/summary_plots_drafts/hits_TF_para_pairs_p_fc.csv', quote = F, row.names = F)
write.csv(nonhits_TF_para_pairs_anypara, file = 'de_analysis/summary_plots_drafts/nonhits_TF_para_pairs_anypara.csv', quote = F, row.names = F)

hits_TF_para_pairs_regulons <- as_tibble(read.csv('supp_analyses/shared_regulons/bulk/hits_tf_analysis.csv', stringsAsFactors = F)) %>%
  inner_join(hits_TF_para_pairs_p_fc %>% dplyr::select(KO_Gene, GEO_ID) %>% dplyr::rename(source.x = KO_Gene))
nonhits_TF_para_pairs_regulons <- as_tibble(read.csv('supp_analyses/shared_regulons/bulk/non-hits_tf_analysis.csv', stringsAsFactors = F)) %>%
  inner_join(nonhits_TF_para_pairs_anypara %>% dplyr::select(KO_Gene, GEO_ID) %>% dplyr::rename(source.x = KO_Gene))


hit_reg_DESeq <- list()
for(ds in unique(hits_TF_para_pairs_regulons$GEO_ID)) {
  
  tdeseq <- as_tibble(read.csv(paste0('deseq_files/', ds, '/differentialExpression_DESeq_allTargets.csv'), stringsAsFactors = F))
  tdeseq$baseMean <- as.numeric(tdeseq$baseMean)
  tdeseq$log2FoldChange <- as.numeric(tdeseq$log2FoldChange)
  tdeseq$lfcSE <- as.numeric(tdeseq$lfcSE)
  tdeseq$stat <- as.numeric(tdeseq$stat)
  tdeseq$pvalue <- as.numeric(tdeseq$pvalue)
  tdeseq$padj <- as.numeric(tdeseq$padj)
  
  temppairs <- hits_TF_para_pairs_regulons %>% 
    filter(GEO_ID == ds) %>%
    dplyr::select(source.x, source.y) %>%
    unique()
  temppairs$tid <- 1:nrow(temppairs)
  
  for(kop in 1:nrow(temppairs)) {
    
    ko = temppairs[kop, "source.x"]
    pa = temppairs[kop, "source.y"]
    
    temp_reg <- hits_TF_para_pairs_regulons %>%
      filter(source.x == ko$source.x,
             source.y == pa$source.y,
             GEO_ID == ds)
    
    
    reg_DEres <- tdeseq %>%
      filter(sampleKO == ko$source.x,
             gene_name %in% temp_reg$target) %>%
      mutate(paralog = pa$source.y,
             GEO_ID = ds)
    
    if(is.null(dim(hit_reg_DESeq))) {
      hit_reg_DESeq <- reg_DEres
    } else {
      hit_reg_DESeq %<>% bind_rows(reg_DEres)
    }
    
  }
}

nonhit_reg_DESeq <- list()
for(ds in unique(nonhits_TF_para_pairs_regulons$GEO_ID)) {
  
  tdeseq <- as_tibble(read.csv(paste0('deseq_files/', ds, '/differentialExpression_DESeq_allTargets.csv'), stringsAsFactors = F))
  tdeseq$baseMean <- as.numeric(tdeseq$baseMean)
  tdeseq$log2FoldChange <- as.numeric(tdeseq$log2FoldChange)
  tdeseq$lfcSE <- as.numeric(tdeseq$lfcSE)
  tdeseq$stat <- as.numeric(tdeseq$stat)
  tdeseq$pvalue <- as.numeric(tdeseq$pvalue)
  tdeseq$padj <- as.numeric(tdeseq$padj)
  
  temppairs <- nonhits_TF_para_pairs_regulons %>% 
    filter(GEO_ID == ds) %>%
    dplyr::select(source.x, source.y) %>%
    unique()
  temppairs$tid <- 1:nrow(temppairs)
  
  for(kop in 1:nrow(temppairs)) {
    
    ko = temppairs[kop, "source.x"]
    pa = temppairs[kop, "source.y"]
    
    temp_reg <- nonhits_TF_para_pairs_regulons %>%
      filter(source.x == ko$source.x,
             source.y == pa$source.y,
             GEO_ID == ds)
    
    
    reg_DEres <- tdeseq %>%
      filter(sampleKO == ko$source.x,
             gene_name %in% temp_reg$target) %>%
      mutate(paralog = pa$source.y,
             GEO_ID = ds)
    
    if(is.null(dim(nonhit_reg_DESeq))) {
      nonhit_reg_DESeq <- reg_DEres
    } else {
      reg_DEres %<>%
        filter(!(paste0(gene_name,GEO_ID) %in% paste0(nonhit_reg_DESeq$gene_name, nonhit_reg_DESeq$GEO_ID))) # remove duplicate regulon members already checked in same dataset
      nonhit_reg_DESeq %<>% bind_rows(reg_DEres)
    }
    
  }
}
regplotdat <- bind_rows(hit_reg_DESeq %>% mutate(isHit='with-TA'),
                        nonhit_reg_DESeq %>% mutate(isHit = 'no-TA'))
set.seed(7492)
bulk_regulon_swarmbox<-ggplot() + 
  geom_beeswarm(data = regplotdat,
                aes(isHit, log2FoldChange), stroke = 0, alpha = 0.3) +
  geom_boxplot(data = regplotdat,
                aes(isHit, log2FoldChange), alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_hline(yintercept = -0.5, linetype = 2) +
  theme_classic() +
  xlab('Upstream TF Type') +
  ylab('Downstream gene log2(fold change)') +
  ggtitle('Regulon expression changes after TF knockout\nBulk RNA-seq CRISPR targets')
ggsave(bulk_regulon_swarmbox, file = 'de_analysis/summary_plots_drafts/bulk_regulon_swarmbox.svg', width = 3, height = 4)

fracDE_hit<-sum(hit_reg_DESeq$padj<0.05 & abs(hit_reg_DESeq$log2FoldChange) > 0.5 , na.rm = T)/length(hit_reg_DESeq$padj)  
fracDE_nonhit<-sum(nonhit_reg_DESeq$padj<0.05 & abs(nonhit_reg_DESeq$log2FoldChange) > 0.5, na.rm = T)/length(nonhit_reg_DESeq$padj)  

frac_DE_tbl <- tibble(
  upstream_TF = c('hit', 'non-hit'),
  frac_regulon_sigDE = c(fracDE_hit, fracDE_nonhit)
) 
write.csv(frac_DE_tbl, file = 'de_analysis/summary_plots_drafts/bulk_fracDE_tbl.csv', row.names = F, quote = F)

FET_dat <- regplotdat %>%
  mutate(isSigDE = !is.na(log2FoldChange) & !is.na(padj) & abs(log2FoldChange)>0.5 & padj<0.05) %>%
  group_by(isHit,isSigDE) %>%
  summarise(nGenes = length(isHit)) %>%
  pivot_wider(names_from = isSigDE, values_from=nGenes)

FET_dat <- as.data.frame(FET_dat)
row.names(FET_dat) <- FET_dat$isHit
FET_dat$isHit<-NULL
fet<-fisher.test(FET_dat)

fet_tbl <- tibble(
  fet_pval = fet$p.value
)
write.csv(fet_tbl, file = 'de_analysis/summary_plots_drafts/bulk_FET_tbl.csv', row.names = F, quote = F)

cv_test_result_fc <- with(regplotdat %>% filter(!is.na(log2FoldChange)), asymptotic_test(2^log2FoldChange, isHit))
cv_test_result_lfc <- with(regplotdat %>% filter(!is.na(log2FoldChange)), asymptotic_test(log2FoldChange, isHit))

cvtesttbl <- tibble(
  xform = c('lfc', 'fc-recommended'),
  pval = c(cv_test_result_lfc$p_value, cv_test_result_fc$p_value)
)
write.csv(cvtesttbl, file = 'de_analysis/summary_plots_drafts/bulk_CVtest_tbl.csv', row.names = F, quote = F)

