library(tidyverse)
library(readxl)
library(gprofiler2)

# genehancer parsed

genehancer_5_19_assoc <- as_tibble(read.table(file = '~/Documents/GeneHancer_AnnotSV_gene_association_scores_v5.19.txt', header = T, stringsAsFactors = F))
genehancer_5_19_elts <- as_tibble(read.table(file = '~/Documents/GeneHancer_AnnotSV_elements_v5.19.txt', header = T, stringsAsFactors = F))

# load CRISPR target genes

bulk_all <- as_tibble(read.csv('/Users/ianmellis/Library/CloudStorage/GoogleDrive-ian.mellis@gmail.com/.shortcut-targets-by-id/1gc1N1Vtn5NLTLibLgUCbN2XpywoLTTgM/MellisI_GoyalY_Collaboration/submissions/submission2_Presub/GenomeBiology/TableS2.csv', header = T, stringsAsFactors = F)) 
bulk_hits <- bulk_all %>%
  filter(bootstrap_p.val < 0.1)
sc_top100s <- as_tibble(read.csv('/Users/ianmellis/Library/CloudStorage/GoogleDrive-ian.mellis@gmail.com/.shortcut-targets-by-id/1gc1N1Vtn5NLTLibLgUCbN2XpywoLTTgM/MellisI_GoyalY_Collaboration/submissions/submission2_Presub/GenomeBiology/TableS3.csv', stringsAsFactors = F, header = T))
sc_all <- 
hits <- c(bulk_hits$KO_Gene, unique(sc_top100s$CRISPR_target_GS))

# load paralogs
#secondary conditions (not for main analysis)


alt_sets <- c('GSE130969-2', 'GSE92872-2', 'GSE145653-2', 'GSE161466-2', 'GSE161466-3', 'GSE161466-4', 'GSE175787-3', 'GSE175787-4', 'GSE175787-5', 'GSE175787-6')

grn_nitc_path <- "~/code/GitHub/grn_nitc/"
dataset_meta <-  as_tibble(read.csv(paste0(grn_nitc_path, 'rnaseq/annotations/dataset_metadata.csv'), stringsAsFactors = F, header = T))
allData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path, 'rnaseq/supp_analyses/length/data/allData.csv')), stringsAsFactors = F, header = T))
filteredData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/filteredData.csv')), stringsAsFactors = F, header = T))
mainData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/mainData.csv')), stringsAsFactors = F, header = T)) %>%
  filter(!(dataset %in% alt_sets)) %>% dplyr::rename(GEO_ID = dataset) %>%
  left_join(dataset_meta, by = 'GEO_ID') %>%
  filter(species == 'human')

#convert to GeneCards IDs for GeneHancer
conv_GCIDs_targs <- gconvert(
  query = mainData$Gene,
  organism = "hsapiens",
  target = "GENECARDS",
  numeric_ns = "",
  mthreshold = 1,
  filter_na = FALSE
)
conv_GCIDs_paras <- gconvert(
  query = mainData$Paralog,
  organism = "hsapiens",
  target = "GENECARDS",
  numeric_ns = "",
  mthreshold = 1,
  filter_na = FALSE
)

mainData$target_GCID <- conv_GCIDs_targs$target
mainData$paralog_GCID <- conv_GCIDs_paras$target

#p65 lookup failed. All paralogs successfully found. Manually found GCID for p65 online: ENSG00000233645
mainData %<>%
  mutate(target_GCID = ifelse(Gene == 'p65', 'ENSG00000233645', target_GCID))
targets <- unique(mainData$target_GCID)


# search for enhancers and promoters associated with CRISPR targets
targ_para_shares_all <- list()
gh_sub_targ_total <- list()
for (gs in targets) {
  
  gh_sub_targ <- genehancer_5_19_assoc %>%
    filter(symbol == gs)
  
  para_sub_targ <- mainData %>%
    filter(target_GCID == gs)
  
  gh_sub_targ_all <- genehancer_5_19_assoc %>%
    filter(GHid %in% gh_sub_targ$GHid) %>%
    inner_join(genehancer_5_19_elts, by = 'GHid') %>%
    mutate(target_GCID = gs) %>%
    dplyr::rename(paralog_GCID = symbol)
  
  if(is.null(gh_sub_targ_total)) {
    gh_sub_targ_total <- gh_sub_targ_all
  } else {
    gh_sub_targ_total %<>% bind_rows(gh_sub_targ_all)
  }
  
  # tmp_paralogs_in_gh_sub_enh <- para_sub_targ %>%
  #   mutate(symbol = paralog_GCID) %>%
  #   inner_join(gh_sub_targ_all, by = 'symbol') %>%
  #   filter(regulatory_element_type == 'Enhancer')
  # 
  # tmp_paralogs_in_gh_sub_pro <- para_sub_targ %>%
  #   mutate(symbol = paralog_GCID) %>%
  #   inner_join(gh_sub_targ_all, by = 'symbol') %>%
  #   filter(regulatory_element_type == 'Promoter')
  # 
}

Allparalogs_summaryshare <- mainData %>%
  left_join(gh_sub_targ_total, by = c('target_GCID', 'paralog_GCID')) %>%
  mutate(shares_enh = !is.na(GHid)) %>%
  dplyr::select(target_GCID, paralog_GCID, FC2, padj, shares_enh) %>%
  group_by(target_GCID, paralog_GCID) %>%
  filter(!is.na(padj)) %>%
  summarise(shares_an_enh = sum(shares_enh) > 0, is_upreg = sum(FC2>0.5 & padj < 0.05) > 0) %>% 
  group_by(shares_an_enh, is_upreg) %>% 
  summarise(n_paralogs = n()) %>%
  pivot_wider(names_from = shares_an_enh, values_from = n_paralogs) %>%
  dplyr::rename(shares_elt = `TRUE`,
                doesnt_share_elt = `FALSE`)

Paralog_summaryshare_mat <- as.matrix(Allparalogs_summaryshare %>% dplyr::select(-is_upreg))

paralog_fet<-fisher.test(Paralog_summaryshare_mat)


# consider all genes, not just paralogs. Lots more data, so will need to iterate through targets and load respective datasets
# ignore gene names without gprofiler-found GeneCards symbols

dataset_meta

unique_sets <- unique(mainData$GEO_ID)

targ_anygene_shares_all <- list()
gh_sub_targ_total_anygene <- list()
anygenes_summaryshare <- list()
for (ds in unique_sets) {
  
  cat('Loading', ds, 'and converting gene IDs...')
  allDES_tmp <- read.csv(paste0(grn_nitc_path, '/rnaseq/deseq_files/', ds, '/differentialExpression_DESeq_allTargets.csv'), stringsAsFactors = F, header = T)
  
  tmp_main <- mainData %>%
    filter(GEO_ID == ds)
  
  unique_targs_in_ds <- unique(tmp_main$target_GCID)
  
  if (ds != 'GSE151825') {
    conv_GCIDs_anygene_tmp <- gconvert(
      query = allDES_tmp$gene_name,
      organism = "hsapiens",
      target = "GENECARDS",
      numeric_ns = "",
      mthreshold = 1,
      filter_na = FALSE
    )
    
    allDES_tmp$othergene_GCID <- conv_GCIDs_anygene_tmp$target
  }
  cat('Working on', unique_targs_in_ds, '...\n')
  
  for (gs in unique_targs_in_ds) {
    
    if (ds == 'GSE151825') {
      cat('Working on', gs, '\n')
      
      allDES_tmp1 <- allDES_tmp %>%
        filter(sampleKO == gs)
      
      conv_GCIDs_anygene_tmp <- gconvert(
        query = allDES_tmp1$gene_name,
        organism = "hsapiens",
        target = "GENECARDS",
        numeric_ns = "",
        mthreshold = 1,
        filter_na = FALSE
      )
      
      allDES_tmp1$othergene_GCID <- conv_GCIDs_anygene_tmp$target
    }
    
    conv_tmp <- conv_GCIDs_targs %>%
      filter(target == gs) %>%
      dplyr::select(input, target) %>% unique()
    
    gs_hgnc <- conv_tmp$input[1]
    
    if (ds == "GSE151825") {
      mainData_tmp <- allDES_tmp1 %>%
        filter(sampleKO == gs_hgnc) %>%
        mutate(target_GCID = gs)
    } else {
      mainData_tmp <- allDES_tmp %>%
        filter(sampleKO == gs_hgnc) %>%
        mutate(target_GCID = gs)
    }
    
    gh_sub_targ <- genehancer_5_19_assoc %>%
      filter(symbol == gs)
    
    all_sub_targ <- mainData %>%
      filter(target_GCID == gs)
    
    gh_sub_targ_all <- genehancer_5_19_assoc %>%
      filter(GHid %in% gh_sub_targ$GHid) %>%
      inner_join(genehancer_5_19_elts, by = 'GHid') %>%
      mutate(target_GCID = gs) %>%
      dplyr::rename(othergene_GCID = symbol)
    
    if(is.null(gh_sub_targ_total_anygene)) {
      gh_sub_targ_total_anygene <- gh_sub_targ_all
    } else {
      gh_sub_targ_total_anygene %<>% bind_rows(gh_sub_targ_all)
    }
    
    anygenes_summaryshare_tmp <- as_tibble(mainData_tmp) %>%
      left_join(gh_sub_targ_all, by = c('target_GCID', 'othergene_GCID')) %>%
      mutate(shares_enh = !is.na(GHid)) %>%
      dplyr::select(target_GCID, othergene_GCID, log2FoldChange, padj, shares_enh) %>%
      group_by(target_GCID, othergene_GCID) %>%
      filter(!is.na(padj)) %>%
      summarise(shares_an_enh = sum(shares_enh) > 0, is_upreg = sum(log2FoldChange>0.5 & padj < 0.05) > 0) %>% 
      group_by(shares_an_enh, is_upreg) %>% 
      summarise(n_genes = length(is_upreg))%>%
      pivot_wider(names_from = shares_an_enh, values_from = n_genes, values_fill = 0) %>%
      bind_rows(tibble(is_upreg = logical(),
                      `FALSE` = integer(),
                      `TRUE` = integer())) %>%
      dplyr::rename(shares_elt = `TRUE`,
                    doesnt_share_elt = `FALSE`) %>%
      mutate(shares_elt = ifelse(is.na(shares_elt), 0, shares_elt),
             target_GHID = gs,
             dataset = ds) 
    
    if(is.null(dim(anygenes_summaryshare))) {
      anygenes_summaryshare <- anygenes_summaryshare_tmp
    } else {
      anygenes_summaryshare %<>% bind_rows(anygenes_summaryshare_tmp)
    }
    
  }
}

unique_experiments <- anygenes_summaryshare %>%
  dplyr::select(dataset, target_GHID) %>%
  unique()
anygene_FET_tbl <- list()
for (r in 1:nrow(unique_experiments)) {
  
  tmp_exp <- unique_experiments[r,]
  
  tmp_sumshare <- anygenes_summaryshare %>%
    filter(dataset == tmp_exp$dataset[1],
           target_GHID == tmp_exp$target_GHID[1])
  
  Anygene_summaryshare_mat <- as.matrix(tmp_sumshare %>% dplyr::select(-c(is_upreg, target_GHID, dataset)))
  
  anygene_fet<-fisher.test(Anygene_summaryshare_mat)
  
  tmp_exp %<>% mutate(FET_pval = anygene_fet$p.value)
  
  if(is.null(dim(anygene_FET_tbl))) {
    anygene_FET_tbl <- tmp_exp
  } else {
    anygene_FET_tbl %<>% bind_rows(tmp_exp)
  }
  
}

write.csv(anygene_FET_tbl, file = paste0(grn_nitc_path, 'rnaseq/de_analysis/summary_plots_drafts/genehancer_anygene_FET.csv'), quote = F, row.names = F)

  # check scores TBD
  for (pgs in para_sub_targ$paralog_GCID) {
    
    gh_sub_targ_para_enh <- gh_sub_targ_all %>%
      filter(symbol == pgs, regulatory_element_type == 'Enhancer') %>%
      summarise(CRISPR_target_GS = gs,
                paralog_GS = pgs,
                n_shared_enhancers = length(symbol),
                avg_score_enh = ifelse(n_shared_enhancers > 0, mean(combined_score), 0))
    
    gh_sub_targ_para_pro <- gh_sub_targ_all %>%
      filter(symbol == pgs, regulatory_element_type == 'Promoter') %>%
      summarise(CRISPR_target_GS = gs,
                paralog_GS = pgs,
                n_shared_promoters = length(symbol),
                avg_score_pro = ifelse(n_shared_promoters > 0, mean(combined_score), 0))
    
    targ_para_temp <- inner_join(gh_sub_targ_para_enh, gh_sub_targ_para_pro, by = c('CRISPR_target_GS', 'paralog_GS'))
    
    if(is.null(dim(targ_para_shares_all))) {
      targ_para_shares_all <- targ_para_temp
    } else {
      targ_para_shares_all %<>% 
        bind_rows(targ_para_shares_all, targ_para_temp)
    }
    
  }
  




# also pull out each other gene