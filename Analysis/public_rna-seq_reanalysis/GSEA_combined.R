## GSEA for hits
library(clusterProfiler)
library(org.Hs.eg.db)

datadir_bulk <- '~/code/grn_nitc/rnaseq/de_analysis/summary_plots_drafts/'
datadir_singlecell <- '/Volumes/IAMYG2/grn_nitc_data/Frangieh2021/stat_graphs/'

#load bulk table
bulk_hits <- as_tibble(read.csv(paste0(datadir_bulk, 'bulk_targets_hits.csv'), header = T, stringsAsFactors = F))
bulk_alltargs <- as_tibble(read.csv(paste0(datadir_bulk, 'bulk_alltargets.csv'), header = T, stringsAsFactors = F))

#consider only human genes
bulk_hits_human <- bulk_hits %>%
  filter(species == 'human')

bulk_hits_mouse <- bulk_hits %>%
  filter(species == 'mouse')

bulk_alltargs_human <- bulk_alltargs %>%
  filter(species == 'human')

bulk_alltargs_mouse <- bulk_alltargs %>%
  filter(species == 'mouse')

#load single-cell table
sc_hits <- as_tibble(read.csv(paste0(datadir_singlecell, 'singlecell_targets_manycells_hits.csv'), header = T, stringsAsFactors = F))
sc_alltargs <- as_tibble(read.csv(paste0(datadir_singlecell, 'singlecell_alltargets_manycells.csv'), header = T, stringsAsFactors = F))

#combine
all_hits <- bind_rows(bulk_res %>% dplyr::select(CRISPR_target_GS),
                      bulk_alltargs %>% dplyr::select(CRISPR_target_GS))
all_alltargs <- 

#clusterProfiler

bulk_hits_human_BP <- enrichGO(gene = bulk_hits_human$CRISPR_target_GS,
                                universe = bulk_alltargs_human$CRISPR_target_GS,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = "SYMBOL",
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                minGSSize = 2)

bulk_hits_allortho_BP <- enrichGO(gene = all_hits$CRISPR_target_GS,
                                universe = all_alltargs$CRISPR_target_GS,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = "SYMBOL",
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05)

sc_hits_BP <- enrichGO(gene = sc_hits$CRISPR_target_GS,
                               universe = sc_alltargs$CRISPR_target_GS,
                               OrgDb         = org.Hs.eg.db,
                               keyType       = "SYMBOL",
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.1,
                               qvalueCutoff  = 0.1,
                               minGSSize = 2)

bhuman_sc_hits_BP <- enrichGO(gene = c(bulk_hits_human$CRISPR_target_GS, sc_hits$CRISPR_target_GS),
                       universe = c(bulk_alltargs_human$CRISPR_target_GS, sc_alltargs$CRISPR_target_GS),
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "SYMBOL",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1,
                       minGSSize = 2)

bhuman_sc_hits_MF <- enrichGO(gene = c(bulk_hits_human$CRISPR_target_GS, sc_hits$CRISPR_target_GS),
                              universe = c(bulk_alltargs_human$CRISPR_target_GS, sc_alltargs$CRISPR_target_GS),
                              OrgDb         = org.Hs.eg.db,
                              keyType       = "SYMBOL",
                              ont           = "MF",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.1,
                              qvalueCutoff  = 0.1,
                              minGSSize = 2)
