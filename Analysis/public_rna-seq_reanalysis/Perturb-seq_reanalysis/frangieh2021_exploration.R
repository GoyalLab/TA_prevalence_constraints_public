library(tidyverse)
library(biomaRt)
library(magrittr)
library(cvequality)

source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')

datadir <- '/Volumes/IAMYG2/grn_nitc_data/' #change as needed
tflist_f <- '~/code/grn_nitc/Resources/Homo_sapiens_TF.txt'
tflist <- as_tibble(read.table(tflist_f, header = T, sep = '\t')) %>%
  dplyr::rename(CRISPR_target_GS = Symbol)

if(!dir.exists(paste0(datadir, 'Frangieh2021/'))){
  dir.create(paste0(datadir, 'Frangieh2021/'))
}

if(!dir.exists(paste0(datadir, 'Frangieh2021/expression_histograms/'))){
  dir.create(paste0(datadir, 'Frangieh2021/expression_histograms/'))
}

if(!dir.exists(paste0(datadir, 'Frangieh2021/stat_histograms/'))){
  dir.create(paste0(datadir, 'Frangieh2021/stat_histograms/'))
}
fra_meta <- as_tibble(read.csv(paste0(datadir, 'RNA_metadata.csv'), stringsAsFactors = F))
fra_meta <- fra_meta[-1,]

#only MOI 1 cells in "Control" condition (no Ifng or coculture), 2000 UMI minimum
fra_meta_moi1_ctl <- fra_meta %>%
  filter(MOI == 1, condition == 'Control', as.numeric(UMI_count) > 2000) 

# as controls use only non-targeting controls, ignoring intergenic-targeting controls
fra_meta_moi1_ctl_NTC <- fra_meta_moi1_ctl %>%
  filter(grepl('NO_SITE', sgRNA))

# get gene-level target for each non-control guide
fra_meta_moi1_ctl_celltargets <- fra_meta_moi1_ctl %>%
  filter(!grepl('NO_SITE', sgRNA), !grepl('GENE_SITE', sgRNA)) %>%
  separate(col = 'sgRNA', into = c('CRISPR_target', 'guideID'), sep='_') %>%
  dplyr::select(NAME, CRISPR_target)

fra_meta_moi1_ctl_targets <- fra_meta_moi1_ctl %>%
  filter(!grepl('NO_SITE', sgRNA), !grepl('GENE_SITE', sgRNA)) %>%
  inner_join(fra_meta_moi1_ctl_celltargets)

crispr_targets <- unique(fra_meta_moi1_ctl_targets$CRISPR_target)

# check that all genes are found as gene symbols in Ensembl
# pull list of paralogs from ensembl
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl") # version 105
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105) # version 105, Dec 2021
geneParaList <- getBM(attributes = c("ensembl_gene_id", 
                                     "external_gene_name",
                                     "hsapiens_paralog_ensembl_gene", 
                                     "hsapiens_paralog_associated_gene_name",
                                     'hsapiens_paralog_subtype',
                                     'hsapiens_paralog_orthology_type',
                                     'hsapiens_paralog_perc_id',
                                     'hsapiens_paralog_perc_id_r1'),
                      filters = 'external_gene_name',
                      values = crispr_targets,
                      mart = human)
crispr_targets[!(crispr_targets %in% geneParaList$external_gene_name)]
# 4 targets not in original guide names: "NUP50-AS1"   "LRRC75A-AS1" "TMEM173"     "ATP5MD" 
# TMEM173 genesymbol is "STING1"
# ATP5MD genesymbol is "ATP5MK"
# NUP50-AS1 genesymbol is "NUP50-DT"
# LRRC75A-AS1 genesymbol is "SNHG29"

# update lists to account for actual gene symbols
fra_meta_moi1_ctl_targets %<>% ungroup() %>%
  mutate(CRISPR_target_GS = case_when(
    CRISPR_target == "TMEM173" ~ 'STING1',
    CRISPR_target == 'ATP5MD' ~ 'ATP5MK',
    CRISPR_target == 'NUP50-AS1' ~ 'NUP50-DT',
    CRISPR_target == 'LRRC75A-AS1' ~ 'SNHG29'
  )) %>%
  mutate(CRISPR_target_GS = ifelse(is.na(CRISPR_target_GS), CRISPR_target, CRISPR_target_GS))

crispr_targets_gs <- unique(fra_meta_moi1_ctl_targets$CRISPR_target_GS)
geneParaList <- getBM(attributes = c("ensembl_gene_id", 
                                     "external_gene_name",
                                     "hsapiens_paralog_ensembl_gene", 
                                     "hsapiens_paralog_associated_gene_name",
                                     'hsapiens_paralog_subtype',
                                     'hsapiens_paralog_orthology_type',
                                     'hsapiens_paralog_perc_id',
                                     'hsapiens_paralog_perc_id_r1'),
                      filters = 'external_gene_name',
                      values = crispr_targets_gs,
                      mart = human)
crispr_targets_gs[!(crispr_targets_gs %in% geneParaList$external_gene_name)] # all accounted for

tempdat_sc <- data.table::fread(paste0(datadir, 'RNA_expression.csv.gz'),header = T, nrows=0, stringsAsFactors = F)

NTCdat_sc <- as_tibble(data.table::fread(paste0(datadir, 'RNA_expression.csv.gz'),
                                          header = T, 
                                          select = c('GENE', fra_meta_moi1_ctl_NTC$NAME), 
                                          stringsAsFactors = F)) 
cat('Data loaded. ')

# for KO cells
stats_sc_targs <- list()
stats_sc_NTCs <- list()
for (crispr_targ_gs in crispr_targets_gs[1:12]) {
  
  cat(paste0('Working on ', crispr_targ_gs, '. '))
  # select samples/guides for gene (all guides)
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == crispr_targ_gs)

  # find its paralogs
  temp_paralogs <- geneParaList %>%
    filter(external_gene_name == crispr_targ_gs)
  
  if(nrow(temp_paralogs) <= 1 | nrow(temp_targ_samps) <= 1) {
    cat('Moving on\n') 
    next
  } 
  
  # load data for the cells of interest and filter to just paralogs of the target, reformat to long
  # theoretically, if methods match the table, these values are ln(TPM+1), by UMI
  tempdat_sc <- as_tibble(data.table::fread(paste0(datadir, 'RNA_expression.csv.gz'),
                                            header = T, 
                                            select = c('GENE', temp_targ_samps$NAME), 
                                            stringsAsFactors = F)) %>%
    inner_join(temp_paralogs %>%
                 dplyr::rename(GENE = hsapiens_paralog_associated_gene_name)) %>%
    pivot_longer(cols = 2:(nrow(temp_targ_samps)+1), names_to = 'NAME', values_to = 'abundance')
  cat('Data loaded. ')
  
  tempdat_sc_NTC <- NTCdat_sc %>%
    inner_join(temp_paralogs %>%
                 dplyr::rename(GENE = hsapiens_paralog_associated_gene_name)) %>%
    pivot_longer(cols = 2:(nrow(fra_meta_moi1_ctl_NTC)+1), names_to = 'NAME', values_to = 'abundance')
  
  cat('Making plots. ')
  ttarg_hists <- ggplot() +
    geom_histogram(data = tempdat_sc %>%
                     inner_join(temp_targ_samps, by = 'NAME'), aes(abundance)) +
    facet_grid(sgRNA~GENE, scales = 'free') +
    theme_classic() +
    ggtitle(paste0('scRNA-seq of cells with guides targeting ', crispr_targ_gs,'\nparalogs of ', crispr_targ_gs))
  
  ggsave(ttarg_hists, file = paste0(datadir, 'Frangieh2021/expression_histograms/', crispr_targ_gs, '_targeting_sgRNAs.pdf'), width = length(unique(tempdat_sc$GENE)), height = length(unique(temp_targ_samps$sgRNA)))
  
  tNTC_hists <- ggplot() +
    geom_histogram(data = tempdat_sc_NTC %>%
                     inner_join(fra_meta_moi1_ctl_NTC, by = 'NAME'), aes(abundance)) +
    facet_grid(sgRNA~GENE, scales = 'free') +
    theme_classic() +
    ggtitle(paste0('scRNA-seq of cells with nontemplate guides\nparalogs of ', crispr_targ_gs))
  
  ggsave(tNTC_hists, file = paste0(datadir, 'Frangieh2021/expression_histograms/', crispr_targ_gs, '_nontargeting_sgRNAs.pdf'), width = length(unique(tempdat_sc_NTC$GENE)), height = length(unique(fra_meta_moi1_ctl_NTC$sgRNA)))
  
  cat('Calculating stats.\n')
  # calculate mean and bimodality coefficient for each guide
  tempstats_sc <- tempdat_sc %>%
    inner_join(temp_targ_samps, by = 'NAME') %>%
    group_by(CRISPR_target_GS, sgRNA, GENE) %>%
    summarise(mean_product = mean(abundance),
              median_product = median(abundance),
              perc_pos = sum(abundance > 0),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              skewness = skewness(abundance + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance))
  
  tempstats_sc_NTC <- tempdat_sc_NTC %>%
    inner_join(fra_meta_moi1_ctl_NTC, by = 'NAME') %>%
    group_by(sgRNA, GENE) %>%
    summarise(mean_product = mean(abundance),
              median_product = median(abundance),
              perc_pos = sum(abundance > 0),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              skewness = skewness(abundance + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance)) %>%
    mutate(CRISPR_target_GS = crispr_targ_gs)
  
  if(is.null(dim(stats_sc_targs))) {
    stats_sc_targs <- tempstats_sc
    stats_sc_NTCs <- tempstats_sc_NTC
  } else {
    stats_sc_targs %<>% bind_rows(tempstats_sc)
    stats_sc_NTCs %<>% bind_rows(tempstats_sc_NTC)
  }
  
}

write.csv(stats_sc_targs, file = paste0(datadir, 'Frangieh2021/stats_sc_targs.csv'), quote = F)
write.csv(stats_sc_NTCs, file = paste0(datadir, 'Frangieh2021/stats_sc_NTCs.csv'), quote = F)

stats_sc_targs_long <-stats_sc_targs %>%
  dplyr::select(CRISPR_target_GS, sgRNA, GENE, mean_product, bimodality_coef) %>%
  pivot_longer(cols = mean_product:bimodality_coef, names_to = 'statistic', values_to = 'value')
  
stats_sc_NTCs_long <- stats_sc_NTCs %>%
  dplyr::select(CRISPR_target_GS, sgRNA, GENE, mean_product, bimodality_coef) %>%
  pivot_longer(cols = mean_product:bimodality_coef, names_to = 'statistic', values_to = 'value')

for (crispr_targ_gs in crispr_targets_gs[1:12]) {
  cat(paste0('Working on ', crispr_targ_gs, '. '))
  # select samples/guides for gene (all guides)
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == crispr_targ_gs)
  
  # find its paralogs
  temp_paralogs <- geneParaList %>%
    filter(external_gene_name == crispr_targ_gs)
  
  if(nrow(temp_paralogs) <= 1 | nrow(temp_targ_samps) <= 1) {
    cat('Moving on\n') 
    next
  } 
  t_NTC <- stats_sc_NTCs_long %>% filter(CRISPR_target_GS == crispr_targ_gs)
  t_tar <- stats_sc_targs_long %>% filter(CRISPR_target_GS == crispr_targ_gs)
  stat_hists <- ggplot() +
    facet_grid(GENE ~ statistic, scales = 'free') +
    geom_histogram(data = t_NTC, aes(x = value), alpha = 0.4) +
    geom_vline(data = t_tar, aes(xintercept = value)) +
    theme_classic() +
    ggtitle(paste0(crispr_targ_gs, ' paralogs mean and BC\nTargeting sgRNAs (lines) vs.\nnon-targeting sgRNAs (histogram)'))
  
  ggsave(stat_hists, file = paste0(datadir, 'Frangieh2021/stat_histograms/', crispr_targ_gs, '_stat_hists.pdf'), width = 4, height = (length(unique(t_NTC$GENE))/2 + 2))
  
}
# look for high mean and high bimodality paralogs among targeted cells vs NTC cells
# plot histogram of each summary stat in NTC and overlay the stat value of each sgRNA

## load stats from YG
td <- paste0(datadir, 'Frangieh2021/fromYG/')
dfs <- list.files(td)
dfps <- paste0(td, dfs)

stats_sc_NTCsYG <- as_tibble(read.csv(dfps[1])) %>%
  bind_rows(as_tibble(read.csv(dfps[2]))) %>%
  bind_rows(as_tibble(read.csv(dfps[3]))) %>%
  bind_rows(as_tibble(read.csv(dfps[4]))) %>%
  bind_rows(as_tibble(read.csv(dfps[5]))) %>%
  bind_rows(as_tibble(read.csv(dfps[6]))) %>% unique()

stats_sc_targsYG <- as_tibble(read.csv(dfps[7])) %>%
  bind_rows(as_tibble(read.csv(dfps[8]))) %>%
  bind_rows(as_tibble(read.csv(dfps[9]))) %>%
  bind_rows(as_tibble(read.csv(dfps[10]))) %>%
  bind_rows(as_tibble(read.csv(dfps[11]))) %>%
  bind_rows(as_tibble(read.csv(dfps[12]))) %>% unique()

stats_sc_NTCsYG_sums <- stats_sc_NTCsYG %>%
  group_by(CRISPR_target_GS, GENE) %>%
  summarise(mean_mean_NTC = mean(mean_product),
            sem_mean_NTC = sd(mean_product)/sqrt(length(mean_product)),
            mean_median_NTC = mean(median_product),
            sem_median_NTC = sd(median_product)/sqrt(length(median_product)),
            mean_percpos_NTC = mean(perc_pos),
            sem_percpos_NTC = sd(perc_pos)/sqrt(length(perc_pos)))

stats_sc_targsYG_sums <- stats_sc_targsYG %>%
  group_by(CRISPR_target_GS, GENE) %>%
  summarise(mean_mean_targ = mean(mean_product),
            sem_mean_targ = sd(mean_product)/sqrt(length(mean_product)),
            mean_median_targ = mean(median_product),
            sem_median_targ = sd(median_product)/sqrt(length(median_product)),
            mean_percpos_targ = mean(perc_pos),
            sem_percpos_targ = sd(perc_pos)/sqrt(length(perc_pos)))



# pull out the most dramatic CRISPR targets by difference between targ and NTC means
top_mean_delta <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         lfc_mean_mean = log2((mean_mean_targ+0.1)/(mean_mean_NTC+0.1))) %>%
  arrange(-delta_mean_mean) %>%
  head(100) %>%
  dplyr::select(CRISPR_target_GS, GENE, mean_mean_targ, mean_mean_NTC, delta_mean_mean, lfc_mean_mean)

top_mean_lfc <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         lfc_mean_mean = log2((mean_mean_targ+0.1)/(mean_mean_NTC+0.1))) %>%
  arrange(-lfc_mean_mean) %>%
  head(100) %>%
  dplyr::select(CRISPR_target_GS, GENE, mean_mean_targ, mean_mean_NTC, delta_mean_mean, lfc_mean_mean)

top_percpos_delta <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_percpos_mean = mean_percpos_targ - mean_percpos_NTC,
         lfc_percpos_mean = log2((mean_percpos_targ+0.01)/(mean_percpos_NTC+0.01))) %>%
  arrange(-delta_percpos_mean) %>%
  head(100) %>%
  dplyr::select(CRISPR_target_GS, GENE, mean_percpos_targ, mean_percpos_NTC, delta_percpos_mean, lfc_percpos_mean)


top_deltaAndLFC_hits <- inner_join(top_mean_lfc, top_mean_delta %>% 
                                     dplyr::select(CRISPR_target_GS, GENE), 
                                   by = c('CRISPR_target_GS', 'GENE'))
stats_sc_NTCsYG_all <- stats_sc_NTCsYG_sums %>%
  inner_join(stats_sc_targsYG_sums, by = c('CRISPR_target_GS', 'GENE')) %>%
  left_join(top_mean_delta %>%
              dplyr::select(CRISPR_target_GS, GENE) %>%
              mutate(istop100_mean = 'TRUE'), by = c('CRISPR_target_GS', 'GENE')) %>%
  mutate(istop100_mean = !is.na(istop100_mean)) %>%
  left_join(top_percpos_delta %>%
              dplyr::select(CRISPR_target_GS, GENE) %>%
              mutate(istop100_percpos = 'TRUE'), by = c('CRISPR_target_GS', 'GENE')) %>%
  mutate(istop100_percpos = !is.na(istop100_percpos)) 

multi_nParas_top50_mean <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         delta_mean_median = mean_median_targ - mean_median_NTC,
         delta_mean_percpos = mean_percpos_targ - mean_percpos_NTC) %>%
  arrange(-delta_mean_mean) %>%
  head(50) %>%
  group_by(CRISPR_target_GS) %>%
  summarise(nParas_top50 = length(CRISPR_target_GS)) %>%
  arrange(-nParas_top50) %>%
  filter(nParas_top50 > 1)

multi_nParas_top50_median <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         delta_mean_median = mean_median_targ - mean_median_NTC,
         delta_mean_percpos = mean_percpos_targ - mean_percpos_NTC) %>%
  arrange(-delta_mean_median) %>%
  head(50) %>%
  group_by(CRISPR_target_GS) %>%
  summarise(nParas_top50 = length(CRISPR_target_GS)) %>%
  arrange(-nParas_top50) %>%
  filter(nParas_top50 > 1)

multi_nParas_top50_percpos <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         delta_mean_median = mean_median_targ - mean_median_NTC,
         delta_mean_percpos = mean_percpos_targ - mean_percpos_NTC) %>%
  arrange(-delta_mean_percpos) %>%
  head(50) %>%
  group_by(CRISPR_target_GS) %>%
  summarise(nParas_top50 = length(CRISPR_target_GS)) %>%
  arrange(-nParas_top50) %>%
  filter(nParas_top50 > 1)

multi_nParas_top50_percpos_paras <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         delta_mean_median = mean_median_targ - mean_median_NTC,
         delta_mean_percpos = mean_percpos_targ - mean_percpos_NTC) %>%
  arrange(-delta_mean_percpos) %>%
  head(50) %>%
  group_by(CRISPR_target_GS) 

set.seed(7648)
perPara_meanVmean_alltargs <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all, aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all, aes(mean_mean_NTC, 
                                                ymin = mean_mean_targ - sem_mean_targ, 
                                                ymax = mean_mean_targ + sem_mean_targ,
                                                color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all, aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                 xmax = mean_mean_NTC + sem_mean_NTC,
                                                 y = mean_mean_targ,
                                                 color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS %in% multi_nParas_top50_mean$CRISPR_target_GS), 
  #                 aes(mean_mean_NTC, mean_mean_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means')

set.seed(7648)
perPara_meanVmean_alltargs_onlytop100 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE'), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_text_repel(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE'), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean, label = GENE), alpha = 0.5) +
  geom_errorbar(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE'), aes(mean_mean_NTC, 
                                                ymin = mean_mean_targ - sem_mean_targ, 
                                                ymax = mean_mean_targ + sem_mean_targ,
                                                color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE'), aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                 xmax = mean_mean_NTC + sem_mean_NTC,
                                                 y = mean_mean_targ,
                                                 color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  theme_classic() +
  ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means')
  

set.seed(7648)
perPara_meanVmean_targsNTCpercpos075 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all %>%
               filter(mean_percpos_NTC > 0.75), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all %>%
                  filter(mean_percpos_NTC > 0.75), aes(mean_mean_NTC, 
                                                ymin = mean_mean_targ - sem_mean_targ, 
                                                ymax = mean_mean_targ + sem_mean_targ,
                                                color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all %>%
                   filter(mean_percpos_NTC > 0.75), aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                 xmax = mean_mean_NTC + sem_mean_NTC,
                                                 y = mean_mean_targ,
                                                 color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS %in% multi_nParas_top50_mean$CRISPR_target_GS), 
  #                 aes(mean_mean_NTC, mean_mean_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells\nPercent-positive 0.75 or higher in controls') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means')

set.seed(7648)
perPara_meanVmean_targsNTCpercpos075_onlytop100 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_text_repel(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean, label = GENE), alpha = 0.5) +
  geom_errorbar(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(mean_mean_NTC, 
                                                                                    ymin = mean_mean_targ - sem_mean_targ, 
                                                                                    ymax = mean_mean_targ + sem_mean_targ,
                                                                                    color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                                                     xmax = mean_mean_NTC + sem_mean_NTC,
                                                                                     y = mean_mean_targ,
                                                                                     color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  theme_classic() +
  ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means')


perPara_medianVmedian_alltargs <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all, aes(mean_median_NTC, mean_median_targ), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all, aes(mean_median_NTC, 
                                                ymin = mean_median_targ - sem_median_targ, 
                                                ymax = mean_median_targ + sem_median_targ), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all, aes(xmin = mean_median_NTC - sem_median_NTC, 
                                                 xmax = mean_median_NTC + sem_median_NTC,
                                                 y = mean_median_targ), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_classic() +
  ggtitle('Per paralog median in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA medians') +
  ylab('Mean of reference-targeting sgRNA medians')

perPara_percposVpercpos_alltargs <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all, aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all, aes(mean_percpos_NTC, 
                                                ymin = mean_percpos_targ - sem_percpos_targ, 
                                                ymax = mean_percpos_targ + sem_percpos_targ, 
                                                color = istop100_percpos), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all, aes(xmin = mean_percpos_NTC - sem_percpos_NTC, 
                                                 xmax = mean_percpos_NTC + sem_percpos_NTC,
                                                 y = mean_percpos_targ, 
                                                 color = istop100_percpos), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
    #                   filter(CRISPR_target_GS == 'TUBB',
    #                          GENE %in% multi_nParas_top50_percpos_paras$GENE), 
    #               aes(mean_percpos_NTC, mean_percpos_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  ggtitle('Per paralog percent-positive in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA percent-positives') +
  ylab('Mean of reference-targeting sgRNA percent-positives') 

perPara_percposVpercpos_alltargs_TUBB <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all, aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all, aes(mean_percpos_NTC, 
                                                ymin = mean_percpos_targ - sem_percpos_targ, 
                                                ymax = mean_percpos_targ + sem_percpos_targ, 
                                                color = istop100_percpos), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all, aes(xmin = mean_percpos_NTC - sem_percpos_NTC, 
                                                 xmax = mean_percpos_NTC + sem_percpos_NTC,
                                                 y = mean_percpos_targ, 
                                                 color = istop100_percpos), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  geom_text_repel(data = stats_sc_NTCsYG_all %>%
                    filter(CRISPR_target_GS == 'TUBB',
                           GENE %in% multi_nParas_top50_percpos_paras$GENE),
                aes(mean_percpos_NTC, mean_percpos_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  ggtitle('Per paralog percent-positive in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA percent-positives') +
  ylab('Mean of reference-targeting sgRNA percent-positives') 

perPara_percposVpercpos_alltargs_onlytop100 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all %>% filter(istop100_percpos == 'TRUE'), aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos), alpha = 0.5, stroke = 0) +
  geom_text_repel(data = stats_sc_NTCsYG_all %>% filter(istop100_percpos == 'TRUE'), aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos, label = GENE), alpha = 0.5) +
  geom_errorbar(data = stats_sc_NTCsYG_all %>% filter(istop100_percpos == 'TRUE'), aes(mean_percpos_NTC, 
                                                ymin = mean_percpos_targ - sem_percpos_targ, 
                                                ymax = mean_percpos_targ + sem_percpos_targ, 
                                                color = istop100_percpos), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all %>% filter(istop100_percpos == 'TRUE'), aes(xmin = mean_percpos_NTC - sem_percpos_NTC, 
                                                 xmax = mean_percpos_NTC + sem_percpos_NTC,
                                                 y = mean_percpos_targ, 
                                                 color = istop100_percpos), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS == 'TUBB',
  #                          GENE %in% multi_nParas_top50_percpos_paras$GENE), 
  #               aes(mean_percpos_NTC, mean_percpos_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  ggtitle('Per paralog percent-positive in\nreference-targeted cells vs. in control cells') +
  xlab('Mean of control sgRNA percent-positives') +
  ylab('Mean of reference-targeting sgRNA percent-positives') 

ggsave(perPara_meanVmean_alltargs, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_alltargs.pdf'), width = 6, height = 6)
ggsave(perPara_meanVmean_alltargs_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_alltargs_onlytop100.pdf'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075.pdf'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075_onlytop100.pdf'), width = 10, height = 10)
ggsave(perPara_medianVmedian_alltargs, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_medianVmedian_alltargs.pdf'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_alltargs, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_alltargs.pdf'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_alltargs_TUBB, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_alltargs_TUBB.pdf'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_alltargs_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_alltargs_onlytop100.pdf'), width = 10, height = 10)

ggsave(perPara_meanVmean_alltargs, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_alltargs.svg'), width = 6, height = 6)
ggsave(perPara_meanVmean_alltargs_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_alltargs_onlytop100.svg'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075.svg'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075_onlytop100.svg'), width = 10, height = 10)
ggsave(perPara_medianVmedian_alltargs, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_medianVmedian_alltargs.svg'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_alltargs, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_alltargs.svg'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_alltargs_TUBB, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_alltargs_TUBB.svg'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_alltargs_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_alltargs_onlytop100.svg'), width = 10, height = 10)

## now filter to reasonable confidence: at least 30 cells across all guides and at least 5 cells per guide for all 3 guides

celldat_all <- list()
for (ct in unique(fra_meta_moi1_ctl_targets$CRISPR_target_GS)) {
  
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == ct)
  
  temp_targ_sampsGs <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == ct) %>%
    group_by(CRISPR_target_GS, sgRNA) %>%
    summarise(nCells_perguide = length(sgRNA)) %>%
    pivot_wider(names_from = sgRNA, values_from = nCells_perguide)
  
  if(ncol(temp_targ_sampsGs) == 3) {temp_targ_sampsGs %<>% mutate(sgRNA3 = 0)}
  if(ncol(temp_targ_sampsGs) == 2) {temp_targ_sampsGs %<>% mutate(sgRNA2 = 0, 
                                                                  sgRNA3 = 0)}
  if(ncol(temp_targ_sampsGs) == 1) {temp_targ_sampsGs %<>% mutate(sgRNA1 = 0,
                                                                  sgRNA2 = 0, 
                                                                  sgRNA3 = 0)}
  
  colnames(temp_targ_sampsGs)[2:4] = paste0('sgRNA', 1:3)
  
  td = tibble(
    CRISPR_target_GS = ct,
    ncells_total = nrow(temp_targ_samps),
    nguides_total = length(unique(temp_targ_samps$sgRNA))
  ) %>%
    inner_join(temp_targ_sampsGs)
  
  if(is.null(dim(celldat_all))){
    celldat_all <- td
  } else {
    celldat_all %<>% bind_rows(td)
  }
}

celldat_manycells <- celldat_all %>%
  filter(ncells_total >= 30,
         sgRNA1 >=5, sgRNA2 >=5, sgRNA3 >=5)

top_mean_delta_manycells <- stats_sc_NTCsYG_all %>%
  inner_join(celldat_manycells, by = 'CRISPR_target_GS') %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         lfc_mean_mean = log2((mean_mean_targ+0.1)/(mean_mean_NTC+0.1))) %>%
  arrange(-delta_mean_mean) %>%
  head(100) %>%
  dplyr::select(CRISPR_target_GS, GENE, mean_mean_targ, mean_mean_NTC, delta_mean_mean, lfc_mean_mean)

top_mean_lfc_manycells <- stats_sc_NTCsYG_all %>%
  inner_join(celldat_manycells, by = 'CRISPR_target_GS') %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         lfc_mean_mean = log2((mean_mean_targ+0.1)/(mean_mean_NTC+0.1))) %>%
  arrange(-lfc_mean_mean) %>%
  head(100) %>%
  dplyr::select(CRISPR_target_GS, GENE, mean_mean_targ, mean_mean_NTC, delta_mean_mean, lfc_mean_mean)

top_percpos_delta_manycells <- stats_sc_NTCsYG_all %>%
  inner_join(celldat_manycells, by = 'CRISPR_target_GS') %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_percpos_mean = mean_percpos_targ - mean_percpos_NTC,
         lfc_percpos_mean = log2((mean_percpos_targ+0.01)/(mean_percpos_NTC+0.01))) %>%
  arrange(-delta_percpos_mean) %>%
  head(100) %>%
  dplyr::select(CRISPR_target_GS, GENE, mean_percpos_targ, mean_percpos_NTC, delta_percpos_mean, lfc_percpos_mean)

stats_sc_NTCsYG_all_manycells <- stats_sc_NTCsYG_sums %>%
  inner_join(stats_sc_targsYG_sums, by = c('CRISPR_target_GS', 'GENE')) %>%
  inner_join(celldat_manycells, by = 'CRISPR_target_GS') %>%
  left_join(top_mean_delta_manycells %>%
              dplyr::select(CRISPR_target_GS, GENE, delta_mean_mean) %>%
              mutate(istop100_mean = 'TRUE'), by = c('CRISPR_target_GS', 'GENE')) %>%
  mutate(istop100_mean = !is.na(istop100_mean)) %>%
  left_join(top_percpos_delta_manycells %>%
              dplyr::select(CRISPR_target_GS, GENE, delta_percpos_mean) %>%
              mutate(istop100_percpos = 'TRUE'), by = c('CRISPR_target_GS', 'GENE')) %>%
  mutate(istop100_percpos = !is.na(istop100_percpos)) 

multi_top100_bymean_manycells_manycells <- stats_sc_NTCsYG_all_manycells %>%
  inner_join(celldat_manycells, by = 'CRISPR_target_GS') %>%
  filter(istop100_mean) %>%
  dplyr::select(CRISPR_target_GS) %>%
  summarise(nt100 = length(CRISPR_target_GS)) %>%
  arrange(-nt100) %>%
  filter(nt100 > 1)
  
set.seed(7648)
perPara_percposVpercpos_manycells <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all_manycells, aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all_manycells, aes(mean_percpos_NTC, 
                                                ymin = mean_percpos_targ - sem_percpos_targ, 
                                                ymax = mean_percpos_targ + sem_percpos_targ, 
                                                color = istop100_percpos), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all_manycells, aes(xmin = mean_percpos_NTC - sem_percpos_NTC, 
                                                 xmax = mean_percpos_NTC + sem_percpos_NTC,
                                                 y = mean_percpos_targ, 
                                                 color = istop100_percpos), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS == 'TUBB',
  #                          GENE %in% multi_nParas_top50_percpos_paras$GENE), 
  #               aes(mean_percpos_NTC, mean_percpos_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  # ggtitle('Per paralog percent-positive in\nreference-targeted cells vs. in control cells\nonly targets with many cells') +
  xlab('Mean of control sgRNA percent-positives') +
  ylab('Mean of reference-targeting sgRNA percent-positives')  +
  theme(legend.position = 'none')

set.seed(7648)
perPara_meanVmean_targsNTCpercpos075_manycells <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all_manycells %>%
               filter(mean_percpos_NTC > 0.75), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all_manycells %>%
                  filter(mean_percpos_NTC > 0.75), aes(mean_mean_NTC, 
                                                       ymin = mean_mean_targ - sem_mean_targ, 
                                                       ymax = mean_mean_targ + sem_mean_targ,
                                                       color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all_manycells %>%
                   filter(mean_percpos_NTC > 0.75), aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                        xmax = mean_mean_NTC + sem_mean_NTC,
                                                        y = mean_mean_targ,
                                                        color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS %in% multi_nParas_top50_mean$CRISPR_target_GS), 
  #                 aes(mean_mean_NTC, mean_mean_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  # ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells\nPercent-positive 0.75 or higher in controls\nonly targets with many cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means') +
  theme(legend.position = 'none')

set.seed(7648)
perPara_meanVmean_manycells <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all_manycells, aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_errorbar(data = stats_sc_NTCsYG_all_manycells, aes(mean_mean_NTC, 
                                                ymin = mean_mean_targ - sem_mean_targ, 
                                                ymax = mean_mean_targ + sem_mean_targ,
                                                color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all_manycells, aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                 xmax = mean_mean_NTC + sem_mean_NTC,
                                                 y = mean_mean_targ,
                                                 color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS %in% multi_nParas_top50_mean$CRISPR_target_GS), 
  #                 aes(mean_mean_NTC, mean_mean_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  # ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells\nonly targets with many cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means') +
  theme(legend.position = 'none')

set.seed(7648)
perPara_meanVmean_manycells_onlytop100 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE'), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_text_repel(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE'), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean, label = GENE), alpha = 0.5) +
  geom_errorbar(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE'), aes(mean_mean_NTC, 
                                                                                    ymin = mean_mean_targ - sem_mean_targ, 
                                                                                    ymax = mean_mean_targ + sem_mean_targ,
                                                                                    color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE'), aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                                                     xmax = mean_mean_NTC + sem_mean_NTC,
                                                                                     y = mean_mean_targ,
                                                                                     color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  theme_classic() +
  # ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells\nonly targets with many cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means') +
  theme(legend.position = 'none')


set.seed(7648)
perPara_meanVmean_targsNTCpercpos075_manycells_onlytop100 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean), alpha = 0.5, stroke = 0) +
  geom_text_repel(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(mean_mean_NTC, mean_mean_targ, color = istop100_mean, label = GENE), alpha = 0.5) +
  geom_errorbar(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(mean_mean_NTC, 
                                                                                                            ymin = mean_mean_targ - sem_mean_targ, 
                                                                                                            ymax = mean_mean_targ + sem_mean_targ,
                                                                                                            color = istop100_mean), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_mean == 'TRUE',mean_percpos_NTC > 0.75), aes(xmin = mean_mean_NTC - sem_mean_NTC, 
                                                                                                             xmax = mean_mean_NTC + sem_mean_NTC,
                                                                                                             y = mean_mean_targ,
                                                                                                             color = istop100_mean), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  theme_classic() +
  # ggtitle('Per paralog mean in\nreference-targeted cells vs. in control cells\nonly targets with many cells') +
  xlab('Mean of control sgRNA means') +
  ylab('Mean of reference-targeting sgRNA means') +
  theme(legend.position = 'none')

perPara_percposVpercpos_manycells_onlytop100 <- ggplot() +
  geom_point(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_percpos == 'TRUE'), aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos), alpha = 0.5, stroke = 0) +
  geom_text_repel(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_percpos == 'TRUE'), aes(mean_percpos_NTC, mean_percpos_targ, color = istop100_percpos, label = GENE), alpha = 0.5) +
  geom_errorbar(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_percpos == 'TRUE'), aes(mean_percpos_NTC, 
                                                                                       ymin = mean_percpos_targ - sem_percpos_targ, 
                                                                                       ymax = mean_percpos_targ + sem_percpos_targ, 
                                                                                       color = istop100_percpos), alpha = 0.5) +
  geom_errorbarh(data = stats_sc_NTCsYG_all_manycells %>% filter(istop100_percpos == 'TRUE'), aes(xmin = mean_percpos_NTC - sem_percpos_NTC, 
                                                                                        xmax = mean_percpos_NTC + sem_percpos_NTC,
                                                                                        y = mean_percpos_targ, 
                                                                                        color = istop100_percpos), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black')) +
  # geom_text_repel(data = stats_sc_NTCsYG_all %>% 
  #                   filter(CRISPR_target_GS == 'TUBB',
  #                          GENE %in% multi_nParas_top50_percpos_paras$GENE), 
  #               aes(mean_percpos_NTC, mean_percpos_targ, label = paste0(GENE,' (', CRISPR_target_GS,')'))) +
  theme_classic() +
  # ggtitle('Per paralog percent-positive in\nreference-targeted cells vs. in control cells\nonly targets with many cells') +
  xlab('Mean of control sgRNA percent-positives') +
  ylab('Mean of reference-targeting sgRNA percent-positives')  +
  theme(legend.position = 'none')

ggsave(perPara_meanVmean_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_manycells.pdf'), width = 6, height = 6)
ggsave(perPara_meanVmean_manycells_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_manycells_onlytop100.pdf'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075_manycells.pdf'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075_manycells_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075_manycells_onlytop100.pdf'), width = 10, height = 10)
ggsave(perPara_percposVpercpos_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_manycells.pdf'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_manycells_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_manycells_onlytop100.pdf'), width = 10, height = 10)

ggsave(perPara_meanVmean_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_manycells.svg'), width = 6, height = 6)
ggsave(perPara_meanVmean_manycells_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_manycells_onlytop100.svg'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075_manycells.svg'), width = 10, height = 10)
ggsave(perPara_meanVmean_targsNTCpercpos075_manycells_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_meanVmean_targsNTCpercpos075_manycells_onlytop100.svg'), width = 10, height = 10)
ggsave(perPara_percposVpercpos_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_manycells.svg'), width = 6, height = 6)
ggsave(perPara_percposVpercpos_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_manycells_wide.svg'), width = 12, height = 6)
ggsave(perPara_percposVpercpos_manycells_onlytop100, file = paste0(datadir, 'Frangieh2021/stat_graphs/perPara_percposVpercpos_manycells_onlytop100.svg'), width = 10, height = 10)


toptargs_manycells <- stats_sc_NTCsYG_all_manycells %>% 
  filter(istop100_percpos == 'TRUE' | (istop100_mean == 'TRUE' & mean_percpos_NTC > 0.75)) %>%
  dplyr::select(CRISPR_target_GS) %>% unique()

write.csv(toptargs_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/singlecell_targets_manycells_hits.csv'), quote = F)
write.csv(stats_sc_NTCsYG_all_manycells %>%
            dplyr::select(CRISPR_target_GS) %>% unique(), file = paste0(datadir, 'Frangieh2021/stat_graphs/singlecell_alltargets_manycells.csv'), quote = F)




topTFs_manycells <- stats_sc_NTCsYG_all_manycells %>% 
  filter(istop100_percpos == 'TRUE' | (istop100_mean == 'TRUE' & mean_percpos_NTC > 0.75)) %>%
  inner_join(tflist %>% dplyr::select(CRISPR_target_GS), by = 'CRISPR_target_GS') %>%
  arrange(-delta_percpos_mean) %>%
  dplyr::select(CRISPR_target_GS, GENE, delta_percpos_mean, delta_mean_mean, mean_percpos_targ, mean_percpos_NTC, mean_mean_targ, mean_mean_NTC) %>%
  dplyr::rename(paralog_gene = GENE)
write.csv(topTFs_manycells, file = paste0(datadir, 'Frangieh2021/stat_graphs/topTFs_manycells.csv'), quote = F)
write.csv(topTFs_manycells %>%
            dplyr::select(CRISPR_target_GS, paralog_gene) %>%
            mutate(downstream_gene = NA,
                   edge_weight_Target = NA,
                   edge_weight_Paralog = NA), file = paste0(datadir, 'Frangieh2021/stat_graphs/topTFs_manycells_forNB.csv'), quote = F)

write.csv(stats_sc_NTCsYG_all_manycells %>% 
            filter(istop100_percpos == 'TRUE' | (istop100_mean == 'TRUE' & mean_percpos_NTC > 0.75)) %>%
            dplyr::rename(Paralog_GS = GENE) %>%
            dplyr::select(CRISPR_target_GS, Paralog_GS, mean_percpos_NTC, mean_percpos_targ, mean_mean_NTC, mean_mean_targ, ncells_total, istop100_percpos, istop100_mean),
          file = '~/code/grn_nitc/fig_drafts_IAM/TableS3.csv')

all_TFs <- stats_sc_NTCsYG_all_manycells %>% 
  inner_join(tflist %>% dplyr::select(CRISPR_target_GS), by = 'CRISPR_target_GS') %>%
  arrange(-delta_percpos_mean) %>%
  dplyr::select(CRISPR_target_GS, GENE, delta_percpos_mean, delta_mean_mean, mean_percpos_targ, mean_percpos_NTC, mean_mean_targ, mean_mean_NTC) %>%
  mutate(delta_mean_mean = ifelse(is.na(delta_mean_mean), mean_mean_targ-mean_mean_NTC, delta_mean_mean),
         delta_percpos_mean = ifelse(is.na(delta_percpos_mean), mean_percpos_targ-mean_percpos_NTC, delta_percpos_mean)) %>%
  dplyr::rename(paralog_gene = GENE)

IQR_l = summary(all_TFs$delta_percpos_mean)['1st Qu.']
IQR_h = summary(all_TFs$delta_percpos_mean)['3rd Qu.']

IQR_nonhitTFs <- all_TFs %>%
  filter(delta_percpos_mean < as.numeric(IQR_h), 
         delta_percpos_mean > as.numeric(IQR_l), 
         mean_percpos_NTC < 0.75,
         !(CRISPR_target_GS %in% topTFs_manycells$CRISPR_target_GS),
         !(paralog_gene %in% topTFs_manycells$CRISPR_target_GS)) 

write.csv(IQR_nonhitTFs %>%
            dplyr::select(CRISPR_target_GS, paralog_gene) %>%
            mutate(downstream_gene = NA,
                   edge_weight_Target = NA,
                   edge_weight_Paralog = NA), file = paste0(datadir, 'Frangieh2021/stat_graphs/nonhitTFs_manycells_forNB.csv'), quote = F)

# pull out the most dramatic
multi_nParas_top50_percpos <- stats_sc_NTCsYG_all %>%
  filter(mean_percpos_targ > 0) %>%
  mutate(delta_mean_mean = mean_mean_targ - mean_mean_NTC,
         delta_mean_median = mean_median_targ - mean_median_NTC,
         delta_mean_percpos = mean_percpos_targ - mean_percpos_NTC) %>%
  arrange(-delta_mean_percpos) %>%
  head(50) %>%
  group_by(CRISPR_target_GS) %>%
  summarise(nParas_top50 = length(CRISPR_target_GS)) %>%
  arrange(-nParas_top50) %>%
  filter(nParas_top50 > 1)

celldat <- list()
for (ct in multi_nParas_top50_percpos$CRISPR_target_GS) {
  
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
        filter(CRISPR_target_GS == ct)
  
  temp_targ_sampsGs <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == ct) %>%
    group_by(CRISPR_target_GS, sgRNA) %>%
    summarise(nCells_perguide = length(sgRNA)) %>%
    pivot_wider(names_from = sgRNA, values_from = nCells_perguide)
  
  if(ncol(temp_targ_sampsGs) == 3) {temp_targ_sampsGs %<>% mutate(sgRNA3 = 0)}
  
  colnames(temp_targ_sampsGs)[2:4] = paste0('sgRNA', 1:3)
  
  td = tibble(
    CRISPR_target_GS = ct,
    ncells_total = nrow(temp_targ_samps),
    nguides_total = length(unique(temp_targ_samps$sgRNA))
  ) %>%
    inner_join(temp_targ_sampsGs)
  
  if(is.null(dim(celldat))){
    celldat <- td
  } else {
    celldat %<>% bind_rows(td)
  }
}

multi_nParas_top50_percpos %>% inner_join(celldat)

## single-cell dists for regulons

# load Nico's table
regulons <- as_tibble(read.csv('~/code/grn_nitc/rnaseq/single-cell/topTFs_manycells_forNB_regulons.csv', stringsAsFactors = F))

NTCdat_sc <- as_tibble(data.table::fread(paste0(datadir, 'RNA_expression.csv.gz'),
                                         header = T, 
                                         select = c('GENE', fra_meta_moi1_ctl_NTC$NAME), 
                                         stringsAsFactors = F)) 
TFdat_sc <- list()
for (targ in unique(regulons$CRISPR_target_GS)) {
  
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == targ)
  
  TFdat_sc[[targ]] <- as_tibble(data.table::fread(paste0(datadir, 'RNA_expression.csv.gz'),
                                                  header = T, 
                                                  select = c('GENE', temp_targ_samps$NAME), 
                                                  stringsAsFactors = F)) 
  
}

TFdat_sc_long <- list()
TFdat_NTC_long <- list()

for (targ in unique(regulons$CRISPR_target_GS)) {
  
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == targ)
  
  gs <- bind_rows(tibble(GENE = targ,
                         type = 'A_CRISPR_target'),
                  regulons %>%
                    filter(CRISPR_target_GS == targ) %>%
                    dplyr::select(paralog_gene) %>%
                    dplyr::rename(GENE = paralog_gene) %>%
                    unique() %>%
                    mutate(type = 'Aprime_Paralog')) %>%
    bind_rows(regulons %>%
                filter(CRISPR_target_GS == targ) %>%
                dplyr::select(downstream_gene) %>%
                dplyr::rename(GENE = downstream_gene) %>%
                unique() %>%
                mutate(type = 'B_Downstream'))
  
  TFdat_sc_long[[targ]] <- TFdat_sc[[targ]] %>%
    inner_join(gs) %>%
    pivot_longer(cols = 2:(nrow(temp_targ_samps)+1), names_to = 'NAME', values_to = 'abundance')
  
  TFdat_NTC_long[[targ]] <- NTCdat_sc %>%
    inner_join(gs) %>%
    pivot_longer(cols = 2:(nrow(fra_meta_moi1_ctl_NTC)+1), names_to = 'NAME', values_to = 'abundance')
  
}

hit_dn_stats <- list()
for (ent in 1:nrow(regulons)) {
  
  tf <- regulons$CRISPR_target_GS[ent]
  pa <- regulons$paralog_gene[ent]
  dn <- regulons$downstream_gene[ent]
  
  temp_targ_vals <- TFdat_sc_long[[tf]] %>%
    filter(GENE %in% c(tf, pa, dn)) %>%
    mutate(cellType = 'targeted')
  
  temp_NTC_vals <- TFdat_NTC_long[[tf]] %>%
    filter(GENE %in% c(tf, pa, dn)) %>%
    mutate(cellType = 'controls')
  
  temp_hists <- ggplot() +
    facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
    geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals), aes(abundance)) +
    theme_classic() +
    geom_text(data = )
  
  ggsave(temp_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/hists_', tf, '_', pa, '_', dn, '.pdf'), width = 6, height = 4)
  ggsave(temp_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/hists_', tf, '_', pa, '_', dn, '.svg'), width = 6, height = 4)
  
  # check change in percpos and mean per downstream target 
  if(nrow(temp_NTC_vals %>% filter(GENE == dn))>0) {
    temp_stats_dn <- temp_NTC_vals %>%
      filter(GENE == dn) %>%
      bind_rows(temp_targ_vals %>% filter(GENE == dn)) %>%
      group_by(GENE, cellType) %>% 
      summarise(mean = mean(abundance),
                percpos = sum(abundance > 0)/length(abundance)) %>%
      pivot_wider(names_from = 'cellType', values_from = c('mean', 'percpos')) %>%
      mutate(delta_mean = mean_targeted - mean_controls,
             delta_percpos = percpos_targeted - percpos_controls)
  
    
    if(is.null(dim(hit_dn_stats))){
      hit_dn_stats = temp_stats_dn
    } else {
      hit_dn_stats %<>% bind_rows(temp_stats_dn)
    }
  }
}
  
tf <- 'SMAD4'
pa <- 'SMAD1'
dn1 <- 'ID3'
dn2 <- 'TNFRSF11B'
dn3 <- 'ID2'

temp_targ_vals <- TFdat_sc_long[[tf]] %>%
  filter(GENE %in% c(tf, pa, dn1, dn2, dn3)) %>%
  mutate(cellType = 'targeted')

temp_NTC_vals <- TFdat_NTC_long[[tf]] %>%
  filter(GENE %in% c(tf, pa, dn1, dn2, dn3)) %>%
  mutate(cellType = 'controls')

sumdat_for5B <- bind_rows(temp_targ_vals, temp_NTC_vals) %>% 
  dplyr::select(GENE, type, cellType) %>%
  unique() %>%
  inner_join(temp_NTC_vals %>%
               filter(GENE %in% c(tf, pa, dn1, dn2, dn3)) %>%
               bind_rows(temp_targ_vals %>% filter(GENE %in% c(tf, pa, dn1, dn2, dn3))) %>%
               group_by(GENE, cellType) %>% 
               summarise(mean = mean(abundance),
                         percpos = sum(abundance > 0)/length(abundance))) 

tfpa_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(tf, pa)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(tf, pa)), aes(x = 2, y = 100, label = paste0('% pos. = ', round(percpos,3))))

ggsave(tfpa_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_tfpa_hists.pdf'), width = 6, height = 4)
ggsave(tfpa_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_tfpa_hists.svg'), width = 6, height = 4)

dn12_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(dn1, dn2)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(dn1, dn2)), aes(x = 2, y = 25, label = paste0('% pos. = ', round(percpos,3))))

ggsave(dn12_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_dn12_hists.pdf'), width = 6, height = 4)
ggsave(dn12_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_dn12_hists.svg'), width = 6, height = 4)

# separate plots for each
tf_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(tf)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(tf)), aes(x = 2, y = 100, label = paste0('% pos. = ', round(percpos,3))))
pa_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(pa)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(pa)), aes(x = 2, y = 100, label = paste0('% pos. = ', round(percpos,3))))
dn1_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(dn1)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(dn1)), aes(x = 2, y = 25, label = paste0('% pos. = ', round(percpos,3))))
dn2_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(dn2)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(dn2)), aes(x = 2, y = 25, label = paste0('% pos. = ', round(percpos,3))))
dn3_hists <- ggplot() +
  facet_grid(cellType~paste0(type,'\n',GENE), scales = 'free_y') +
  geom_histogram(data = bind_rows(temp_targ_vals, temp_NTC_vals) %>% filter(GENE %in% c(dn3)), aes(abundance)) +
  theme_classic() +
  geom_text(data = sumdat_for5B %>% filter(GENE %in% c(dn3)), aes(x = 2, y = 100, label = paste0('% pos. = ', round(percpos,3))))

ggsave(tf_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_tf_hists.svg'), width = 6, height = 4)
ggsave(pa_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_pa_hists.svg'), width = 6, height = 4)
ggsave(dn1_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_dn1_hists.svg'), width = 6, height = 4)
ggsave(dn2_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_dn2_hists.svg'), width = 6, height = 4)
ggsave(dn3_hists, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/SMADs_dn3_hists.svg'), width = 6, height = 4)

ggplot(bind_rows(temp_targ_vals, temp_NTC_vals),
       aes(Aprime_Paralog_SMAD1, B_Downstream_TNFRSF11B, color = A_CRISPR_target_SMAD4==0)) +
  geom_point(position = position_jitter(width = 0.1, height = 0.1, seed = 3476), alpha = 0.5, stroke = 0) + 
  facet_grid(.~cellType) +
  theme_classic()

# write.csv(regulons %>% dplyr::select(-X) %>% dplyr::rename(MOR_TargetTF = edge_weight_Target,
#                                                            MOR_ParalogTF = edge_weight_Paralog), file = '~/code/grn_nitc/fig_drafts_IAM/TableS3.csv')

# load Nico's table - nonhits
regulons_nh <- as_tibble(read.csv('~/code/grn_nitc/rnaseq/single-cell/nonhit_tf_analysis_regulons.csv', stringsAsFactors = F))

TFdat_sc_nh <- list()
for (targ in unique(regulons_nh$CRISPR_target_GS)) {
  
  cat('Working on ', targ, '\n')
  
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == targ)
  
  TFdat_sc_nh[[targ]] <- as_tibble(data.table::fread(paste0(datadir, 'RNA_expression.csv.gz'),
                                                  header = T, 
                                                  select = c('GENE', temp_targ_samps$NAME), 
                                                  stringsAsFactors = F)) 
  
}

TFdat_sc_nh_long <- list()
TFdat_NTC_nh_long <- list()

for (targ in unique(regulons_nh$CRISPR_target_GS)) {
  
  temp_targ_samps <- fra_meta_moi1_ctl_targets %>%
    filter(CRISPR_target_GS == targ)
  
  gs <- bind_rows(tibble(GENE = targ,
                         type = 'A_CRISPR_target'),
                  regulons_nh %>%
                    filter(CRISPR_target_GS == targ) %>%
                    dplyr::select(paralog_gene) %>%
                    dplyr::rename(GENE = paralog_gene) %>%
                    unique() %>%
                    mutate(type = 'Aprime_Paralog')) %>%
    bind_rows(regulons_nh %>%
                filter(CRISPR_target_GS == targ) %>%
                dplyr::select(downstream_gene) %>%
                dplyr::rename(GENE = downstream_gene) %>%
                unique() %>%
                mutate(type = 'B_Downstream'))
  
  TFdat_sc_nh_long[[targ]] <- TFdat_sc_nh[[targ]] %>%
    inner_join(gs) %>%
    pivot_longer(cols = 2:(nrow(temp_targ_samps)+1), names_to = 'NAME', values_to = 'abundance')
  
  TFdat_NTC_nh_long[[targ]] <- NTCdat_sc %>%
    inner_join(gs) %>%
    pivot_longer(cols = 2:(nrow(fra_meta_moi1_ctl_NTC)+1), names_to = 'NAME', values_to = 'abundance')
  
}


# pull out percent pos and means across targeting samples and across NTC samples for each paralog and downstream gene
# check difference. Is the difference minimal? How does the range compare to downstream genes in regulons that aren't hits?


nh_dn_stats <- list()
for (ent in 1:nrow(regulons_nh)) {
  
  tf <- regulons_nh$CRISPR_target_GS[ent]
  pa <- regulons_nh$paralog_gene[ent]
  dn <- regulons_nh$downstream_gene[ent]
  
  temp_targ_vals <- TFdat_sc_nh_long[[tf]] %>%
    filter(GENE %in% c(tf, pa, dn)) %>%
    mutate(cellType = 'targeted')
  
  temp_NTC_vals <- TFdat_NTC_nh_long[[tf]] %>%
    filter(GENE %in% c(tf, pa, dn)) %>%
    mutate(cellType = 'controls')
  
  # check change in percpos and mean per downstream target 
  if(nrow(temp_NTC_vals %>% filter(GENE == dn))>0) {
    temp_stats_dn <- temp_NTC_vals %>%
      filter(GENE == dn) %>%
      bind_rows(temp_targ_vals %>% filter(GENE == dn)) %>%
      group_by(GENE, cellType) %>% 
      summarise(mean = mean(abundance),
                percpos = sum(abundance > 0)/length(abundance)) %>%
      pivot_wider(names_from = 'cellType', values_from = c('mean', 'percpos')) %>%
      mutate(delta_mean = mean_targeted - mean_controls,
             delta_percpos = percpos_targeted - percpos_controls)
    
    
    if(is.null(dim(nh_dn_stats))){
      nh_dn_stats = temp_stats_dn
    } else {
      nh_dn_stats %<>% bind_rows(temp_stats_dn)
    }
  }
}

cv_hit <- hit_dn_stats %>% ungroup() %>%
  mutate(tfType = 'TA')  %>%
  filter(percpos_controls < 0.75) %>%
  summarise(cv_delta_percpos = sd(delta_percpos)/mean(delta_percpos))

cv_nonhit <- nh_dn_stats %>% ungroup() %>%
  mutate(tfType = 'TA')  %>%
  filter(percpos_controls < 0.75) %>%
  summarise(cv_delta_percpos = sd(delta_percpos)/mean(delta_percpos))

dns_delta_percpos_le075 <- bind_rows(hit_dn_stats %>% mutate(tfType = 'TA'),
                                  nh_dn_stats %>% mutate(tfType = 'no-TA') %>%
                                    filter(percpos_controls < 0.75))
set.seed(2189)
dn_delta_percpos_swarmbox <- ggplot() +
  geom_beeswarm(data = dns_delta_percpos_le075,
                aes(x = tfType, y = delta_percpos), alpha = 0.3, stroke = 0) +
  geom_boxplot(data = dns_delta_percpos_le075,
               aes(x = tfType, y = delta_percpos), alpha = 0.5) +
  theme_classic() +
  ylab('Change in percent positive\nDownstream target') +
  xlab('Does upstream TF display\ntranscriptional adaptation')
ggsave(dn_delta_percpos_swarmbox, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/dn_delta_percpos_swarmbox.pdf'), width = 3, height = 4)
ggsave(dn_delta_percpos_swarmbox, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/dn_delta_percpos_swarmbox.svg'), width = 3, height = 4)


# downsample non-hit deltas and recalculate CVs for an empirical CV. 49 values in TA group

set.seed(3489)
nst = 1000
cvs_ds_nh <- tibble(
  index = 1:nst,
  cv_ds = NA,
  sd_ds = NA
)
nh_dn_stats_le075 <- nh_dn_stats %>% filter(percpos_controls < 0.75)
for (ns in 1:nst) {
  
  inds_samp <- sample(1:nrow(nh_dn_stats_le075), 49, replace = F)
  
  temp_samp <- nh_dn_stats_le075[inds_samp,]
  
  cvs_ds_nh$cv_ds[ns] <- sd(temp_samp$delta_percpos)/mean(temp_samp$delta_percpos)
  cvs_ds_nh$sd_ds[ns] <- sd(temp_samp$delta_percpos)
  
}

hit_sds <- hit_dn_stats %>% ungroup() %>%
  mutate(tfType = 'TA')  %>%
  filter(percpos_controls < 0.75) %>%
  summarise(sd_delta_percpos = sd(delta_percpos))
ds_delta_percpos_sds_histline <- ggplot() +
  geom_histogram(data = cvs_ds_nh, aes(sd_ds)) +
  geom_vline(xintercept = hit_sds$sd_delta_percpos, color = 'red') +
  theme_classic() +
  xlab('Standard deviation')
ggsave(ds_delta_percpos_sds_histline, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/ds_delta_percpos_sds_histline.pdf'), width = 4, height = 4)
ggsave(ds_delta_percpos_sds_histline, file = paste0(datadir, 'Frangieh2021/stat_graphs/regulon_graphs/ds_delta_percpos_sds_histline.svg'), width = 4, height = 4)

cv_test_result <- with(bind_rows(hit_dn_stats %>% mutate(tfType = 'TA'),
                              nh_dn_stats %>% mutate(tfType = 'no-TA') %>%
                                    filter(percpos_controls < 0.75)), asymptotic_test(delta_percpos, tfType))
