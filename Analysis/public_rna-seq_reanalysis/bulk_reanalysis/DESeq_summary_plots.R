library(tidyverse)
library(svglite)
library(magrittr)
library(ggrepel)

setwd('~/code/grn_nitc/rnaseq/')

datasets <- as_tibble(read.csv('annotations/dataset_metadata.csv', stringsAsFactors = F, header = T))
datasets2 <- as_tibble(list.files('de_analysis/p_only/percent_upregulated/tabular_format/')) %>%
  separate(value, 'GEO_ID', '.csv')

#secondary conditions (not for main analysis)
alt_sets <- c('GSE130969-2', 'GSE92872-2', 'GSE145653-2', 'GSE161466-2', 'GSE161466-3', 'GSE161466-4', 'GSE175787-3', 'GSE175787-4', 'GSE175787-5', 'GSE175787-6')

datasets3 <- datasets %>% 
  inner_join(datasets2) %>%
  filter(!(GEO_ID %in% alt_sets))

perup_rnaseq_p <- list()
perup_rnaseq_p_fc <- list()
perup_rnaseq_p_fc_b <- list()
fc_rnaseq_p <- list()
fc_rnaseq_p_fc <- list()
fc_rnaseq_p_fc_b <- list()
# for (geo in datasets$GEO_ID) {
for (geo in datasets3$GEO_ID) {
  cat(paste0('Working on ', geo, '\n'))
  temp_perup_p <- as_tibble(read.csv(paste0('de_analysis/p_only/percent_upregulated/tabular_format/', geo,'.csv'), stringsAsFactors = F, header = T)) %>%
    mutate(GEO_ID = geo)
  temp_perup_p_fc <- as_tibble(read.csv(paste0('de_analysis/p_fc/percent_upregulated/tabular_format/', geo,'.csv'), stringsAsFactors = F, header = T)) %>%
    mutate(GEO_ID = geo)
  temp_perup_p_fc_b <- as_tibble(read.csv(paste0('de_analysis/p_fc_basemean/percent_upregulated/tabular_format/', geo,'.csv'), stringsAsFactors = F, header = T)) %>%
    mutate(GEO_ID = geo)
  
  temp_fc_p <- as_tibble(read.csv(paste0('de_analysis/p_only/percent_upregulated/tabular_format/', geo,'.csv'), stringsAsFactors = F, header = T)) %>%
    mutate(GEO_ID = geo)
  temp_fc_p_fc <- as_tibble(read.csv(paste0('de_analysis/p_fc/percent_upregulated/tabular_format/', geo,'.csv'), stringsAsFactors = F, header = T)) %>%
    mutate(GEO_ID = geo)
  temp_fc_p_fc_b <- as_tibble(read.csv(paste0('de_analysis/p_fc_basemean/percent_upregulated/tabular_format/', geo,'.csv'), stringsAsFactors = F, header = T)) %>%
    mutate(GEO_ID = geo)
  
  if (is.null(dim(perup_rnaseq_p))) {
    perup_rnaseq_p <- temp_perup_p
    perup_rnaseq_p_fc <- temp_perup_p_fc
    perup_rnaseq_p_fc_b <- temp_perup_p_fc_b
    
    fc_rnaseq_p <- temp_fc_p
    fc_rnaseq_p_fc <- temp_fc_p_fc
    fc_rnaseq_p_fc_b <- temp_fc_p_fc_b
  } else {
    perup_rnaseq_p %<>% bind_rows(temp_perup_p)
    perup_rnaseq_p_fc %<>% bind_rows(temp_perup_p_fc)
    perup_rnaseq_p_fc_b %<>% bind_rows(temp_perup_p_fc_b)
    
    fc_rnaseq_p %<>% bind_rows(temp_fc_p)
    fc_rnaseq_p_fc %<>% bind_rows(temp_fc_p_fc)
    fc_rnaseq_p_fc_b %<>% bind_rows(temp_fc_p_fc_b)
  }
  
}

perup_summary_p <- perup_rnaseq_p %>%
  dplyr::select(-'X') %>%
  pivot_wider(names_from = variable, values_from = c(Percent.Upregulated, boot_upreg, SE, StdDev, p_value)) %>% 
  dplyr::select(-c("boot_upreg_Paralogs", "SE_Paralogs", "p_value_Paralogs", "boot_upreg_Bootstrapped Genes")) %>%
  mutate(xmax_se = `Percent.Upregulated_Bootstrapped Genes` + `SE_Bootstrapped Genes`,
         xmin_se = `Percent.Upregulated_Bootstrapped Genes` - `SE_Bootstrapped Genes`,
         xmax_sd = `Percent.Upregulated_Bootstrapped Genes` + `StdDev_Bootstrapped Genes`,
         xmin_sd = `Percent.Upregulated_Bootstrapped Genes` - `StdDev_Bootstrapped Genes`) %>%
  mutate(xmin_sd = ifelse(xmin_sd < 0, 0, xmin_sd))
perup_summary_p_fc <- perup_rnaseq_p_fc %>%
  dplyr::select(-'X') %>%
  pivot_wider(names_from = variable, values_from = c(Percent.Upregulated, boot_upreg, SE, StdDev, p_value)) %>% 
  dplyr::select(-c("boot_upreg_Paralogs", "SE_Paralogs", "p_value_Paralogs", "boot_upreg_Bootstrapped Genes")) %>%
  mutate(xmax_se = `Percent.Upregulated_Bootstrapped Genes` + `SE_Bootstrapped Genes`,
         xmin_se = `Percent.Upregulated_Bootstrapped Genes` - `SE_Bootstrapped Genes`,
         xmax_sd = `Percent.Upregulated_Bootstrapped Genes` + `StdDev_Bootstrapped Genes`,
         xmin_sd = `Percent.Upregulated_Bootstrapped Genes` - `StdDev_Bootstrapped Genes`) %>%
  mutate(xmin_sd = ifelse(xmin_sd < 0, 0, xmin_sd))
perup_summary_p_fc_b <- perup_rnaseq_p_fc_b %>%
  dplyr::select(-'X') %>%
  pivot_wider(names_from = variable, values_from = c(Percent.Upregulated, boot_upreg, SE, StdDev, p_value)) %>% 
  dplyr::select(-c("boot_upreg_Paralogs", "SE_Paralogs", "p_value_Paralogs", "boot_upreg_Bootstrapped Genes")) %>%
  mutate(xmax_se = `Percent.Upregulated_Bootstrapped Genes` + `SE_Bootstrapped Genes`,
         xmin_se = `Percent.Upregulated_Bootstrapped Genes` - `SE_Bootstrapped Genes`,
         xmax_sd = `Percent.Upregulated_Bootstrapped Genes` + `StdDev_Bootstrapped Genes`,
         xmin_sd = `Percent.Upregulated_Bootstrapped Genes` - `StdDev_Bootstrapped Genes`) %>%
  mutate(xmin_sd = ifelse(xmin_sd < 0, 0, xmin_sd))

TableS2 <- perup_summary_p_fc %>%
  dplyr::select(KO_Gene, GEO_ID, Percent.Upregulated_Paralogs, `Percent.Upregulated_Bootstrapped Genes`, `p_value_Bootstrapped Genes`) %>%
  dplyr::rename(Observed_fractionUp = Percent.Upregulated_Paralogs,
                Expected_bootstrap_fractionUp = `Percent.Upregulated_Bootstrapped Genes`,
                bootstrap_p.val = `p_value_Bootstrapped Genes`) %>%
  arrange(bootstrap_p.val)
write.csv(TableS2, file = '../fig_drafts_IAM/TableS2.csv')

set.seed(3462)
perup_sum_scatter_p <- ggplot() + 
  geom_point(data = perup_summary_p, aes(`Percent.Upregulated_Bootstrapped Genes`, Percent.Upregulated_Paralogs, color = `p_value_Bootstrapped Genes`)) + 
  geom_errorbarh(data = perup_summary_p, aes(xmin = xmin_sd,
                                             xmax = xmax_sd,
                                             y = Percent.Upregulated_Paralogs,
                                             color = `p_value_Bootstrapped Genes`), height = 0.02) + 
  geom_text_repel(data = perup_summary_p %>% filter(`p_value_Bootstrapped Genes` < 0.1), aes(`Percent.Upregulated_Bootstrapped Genes`, Percent.Upregulated_Paralogs, color = `p_value_Bootstrapped Genes`, label = KO_Gene), seed = 8934) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  scale_color_distiller(palette = 'Spectral', name = 'Bootstrap p.val') +
  xlim(c(0,1)) + 
  xlab('Expected fraction upregulated') +
  ylim(c(0,1)) +
  ylab('Observed fraction upregulated') +
  theme_classic() +
  ggtitle('DESeq2 filters: p.adj < 0.05 only\nLabeled points bootstrap p.val < 0.1')

perup_sum_scatter_p_fc <- ggplot() + 
  geom_point(data = perup_summary_p_fc, aes(`Percent.Upregulated_Bootstrapped Genes`, Percent.Upregulated_Paralogs, color = `p_value_Bootstrapped Genes`)) + 
  geom_errorbarh(data = perup_summary_p_fc, aes(xmin = xmin_sd,
                                                xmax = xmax_sd,
                                                y = Percent.Upregulated_Paralogs,
                                                color = `p_value_Bootstrapped Genes`), height = 0.02) + 
  geom_text_repel(data = perup_summary_p_fc %>% filter(`p_value_Bootstrapped Genes` < 0.1), aes(`Percent.Upregulated_Bootstrapped Genes`, Percent.Upregulated_Paralogs, color = `p_value_Bootstrapped Genes`, label = KO_Gene), seed = 8934) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  scale_color_distiller(palette = 'Spectral', name = 'Bootstrap p.val') +
  xlim(c(0,1)) + 
  xlab('Expected fraction upregulated') +
  ylim(c(-0.03,1.03)) +
  ylab('Observed fraction upregulated') +
  theme_classic() +
  ggtitle('DESeq2 filters: p.adj < 0.05, LFC > 0.5\nLabeled points bootstrap p.val < 0.1')

perup_sum_scatter_p_fc_b <- ggplot() + 
  geom_point(data = perup_summary_p_fc_b, aes(`Percent.Upregulated_Bootstrapped Genes`, Percent.Upregulated_Paralogs, color = `p_value_Bootstrapped Genes`)) + 
  geom_errorbarh(data = perup_summary_p_fc_b, aes(xmin = xmin_sd,
                                                  xmax = xmax_sd,
                                                  y = Percent.Upregulated_Paralogs,
                                                  color = `p_value_Bootstrapped Genes`), height = 0.02) + 
  geom_text_repel(data = perup_summary_p_fc_b %>% filter(`p_value_Bootstrapped Genes` < 0.1), aes(`Percent.Upregulated_Bootstrapped Genes`, Percent.Upregulated_Paralogs, color = `p_value_Bootstrapped Genes`, label = KO_Gene), seed = 8934) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  scale_color_distiller(palette = 'Spectral', name = 'Bootstrap p.val') +
  xlim(c(0,1)) + 
  xlab('Expected fraction upregulated') +
  ylim(c(0,1)) +
  ylab('Observed fraction upregulated') +
  theme_classic() +
  ggtitle('DESeq2 filters: p.adj < 0.05, LFC > 0.5, basemean > 10\nLabeled points bootstrap p.val < 0.1')

ggsave(perup_sum_scatter_p, filename = 'de_analysis/summary_plots_drafts/perup_sum_scatter_p_only.pdf', width = 6, height = 4.5)
ggsave(perup_sum_scatter_p_fc, filename = 'de_analysis/summary_plots_drafts/perup_sum_scatter_p_fc.pdf', width = 6, height = 4.5)
ggsave(perup_sum_scatter_p_fc_b, filename = 'de_analysis/summary_plots_drafts/perup_sum_scatter_p_fc_b.pdf', width = 6, height = 4.5)

ggsave(perup_sum_scatter_p, filename = 'de_analysis/summary_plots_drafts/perup_sum_scatter_p_only.svg', width = 6, height = 4.5)
ggsave(perup_sum_scatter_p_fc, filename = 'de_analysis/summary_plots_drafts/perup_sum_scatter_p_fc.svg', width = 6, height = 4.5)
ggsave(perup_sum_scatter_p_fc_b, filename = 'de_analysis/summary_plots_drafts/perup_sum_scatter_p_fc_b.svg', width = 6, height = 4.5)

qq_perup_p_fc <- ggplot() +
  geom_point(data = perup_summary_p_fc %>% mutate(exp_p = rank(`p_value_Bootstrapped Genes`)/74),
             aes(-log10(exp_p), -log10(`p_value_Bootstrapped Genes`), color = `p_value_Bootstrapped Genes`)) +
  geom_text_repel(data = perup_summary_p_fc %>% mutate(exp_p = rank(`p_value_Bootstrapped Genes`)/74) %>% filter(`p_value_Bootstrapped Genes` < 0.1), 
                  aes(-log10(exp_p), -log10(`p_value_Bootstrapped Genes`), color = `p_value_Bootstrapped Genes`, label = KO_Gene), seed = 8934) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_distiller(palette = 'Spectral', name = 'Bootstrap p.val') +
  theme_classic() +
  ylab('Observed bootstrap -log10(p-value)') +
  xlab('Theoretical -log10(p-value)')

## El-Brolosy targets: Fermt2, Actg1, Actb
elbro <- perup_summary_p_fc %>% 
  filter(KO_Gene %in% c('Fermt2', 'Actg1', 'Actb')) %>%
  pivot_longer(cols = c(Percent.Upregulated_Paralogs, `Percent.Upregulated_Bootstrapped Genes`), names_to = 'statType', values_to = 'meanVal') %>%
  mutate(boot_sd_up = ifelse(statType == 'Percent.Upregulated_Bootstrapped Genes', `StdDev_Bootstrapped Genes`, 0),
         boot_sd_dn = ifelse(statType == 'Percent.Upregulated_Bootstrapped Genes', 
                             ifelse(meanVal<`StdDev_Bootstrapped Genes`, meanVal,`StdDev_Bootstrapped Genes`),
                                    0)) %>%
  dplyr::select(KO_Gene, statType, meanVal, boot_sd_up, boot_sd_dn)

elbro_bars <- ggplot(elbro, aes(x = KO_Gene, weight = meanVal, ymin = meanVal-boot_sd_dn, ymax = meanVal+boot_sd_up, fill = statType)) +
  geom_bar(position = position_dodge(), aes(y=meanVal), stat = 'identity') +
  geom_errorbar(position = position_dodge(width = 0.9), color = 'black', width = 0.25) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab('Fraction of paralogs upregulated') +
  xlab('Knockout target') +
  ggtitle('El-Brolosy targets\nTeal = observed fraction of paralogs upregulated\nRed = expected fraction of paralogs upregulated (bootstrap mean)')
ggsave(elbro_bars,  filename = 'de_analysis/summary_plots_drafts/elbrolosy_targets_bars.svg', width = 6, height = 6)

# isolate hits and present LFCs

hits_p_fc <- perup_summary_p_fc %>% 
  filter(`p_value_Bootstrapped Genes` < 0.1)

write.csv(hits_p_fc %>% dplyr::rename(CRISPR_target_GS = KO_Gene) %>% inner_join(datasets3 %>% dplyr::select(GEO_ID, species), by = 'GEO_ID'), 
          file = 'de_analysis/summary_plots_drafts/bulk_targets_hits.csv',
          quote = F, row.names = F)
write.csv(perup_summary_p_fc %>% dplyr::rename(CRISPR_target_GS = KO_Gene) %>% inner_join(datasets3 %>% dplyr::select(GEO_ID, species), by = 'GEO_ID'), 
          file = 'de_analysis/summary_plots_drafts/bulk_alltargets.csv',
          quote = F, row.names = F)

fc2dat <- list()
corparas <- list()
for (ko in hits_p_fc$KO_Gene) {
  
  t_hit <- hits_p_fc %>% filter(KO_Gene == ko)
  
  tda <- as_tibble(read.csv(paste0('de_analysis/p_fc/FC2/tabular_format/', t_hit$GEO_ID,'.csv'))) %>%
    filter(KO_Gene == ko) %>%
    mutate(GEO_ID = t_hit$GEO_ID) %>%
    dplyr::select(-X) %>%
    filter(!is.na(padj))
  
  if(is.null(dim(fc2dat))) {
    fc2dat <- tda
  } else {
    fc2dat %<>% bind_rows(tda)
  }
  
  if (length(tda$Perc_ID) > 2) {
  cort <- cor.test(tda$Perc_ID, tda$FC2)
  if(is.null(dim(corparas))){
    corparas <- tibble(KO_Gene = ko,
                       Cor = cort$estimate,
                       P.val = cort$p.value)
  } else {
    corparas %<>% bind_rows(tibble(KO_Gene = ko,
                                   Cor = cort$estimate,
                                   P.val = cort$p.value))
  }
  }
}

fc2dat_sum <- fc2dat %>%
  group_by(KO_Gene) %>%
  summarise(mean_FC2 = mean(FC2))

set.seed(6233)
fc2plot_p_fc <- ggplot() +
  geom_jitter(data = fc2dat %>% mutate(sig_DE = padj<0.05), aes(KO_Gene, FC2, color = sig_DE), width = 0.1, height = 0) +
  geom_bar(data = fc2dat_sum, aes(KO_Gene, mean_FC2), alpha = 0.3, stat = 'identity') +
  theme_classic() +
  scale_color_manual(values = c('FALSE' = 'grey50', 'TRUE' = 'black'), name = 'DESeq2 p.adj < 0.05') +
  xlab('CRISPR target') +
  ylab('Paralog log2(fold change)') +
  ggtitle('DESeq2 filters: p.adj < 0.05, LFC > 0.5') +
  geom_hline(yintercept = 0.5, linetype = 2)
ggsave(fc2plot_p_fc, filename = 'de_analysis/summary_plots_drafts/fc2plot_p_fc.pdf', width = 8, height = 5.5)

set.seed(6233)
fc2plot_sDE_p_fc <- ggplot() +
  geom_point(data = fc2dat %>% filter(padj < 0.05, FC2 > 0.5), aes(KO_Gene, FC2)) +
  geom_text_repel(data = fc2dat %>% filter(padj < 0.05, FC2 > 0.5), aes(KO_Gene, FC2, label = Paralog)) +
  theme_classic() +
  xlab('CRISPR target') +
  ylab('Paralog log2(fold change)') +
  ggtitle('DESeq2 filters: p.adj < 0.05, LFC > 0.5') +
  geom_hline(yintercept = 0.5, linetype = 2) +
  ylim(c(0, 10))
ggsave(fc2plot_sDE_p_fc, filename = 'de_analysis/summary_plots_drafts/fc2plot_onlySigDE_p_fc.pdf', width = 8, height = 5.5)


foldchanges_hitsonly_p_fc <- fc2dat %>% filter(padj < 0.05, FC2 > 0.5)
write.csv(foldchanges_hitsonly_p_fc, file ='de_analysis/summary_plots_drafts/foldchanges_hitsonly_p_fc.csv', quote = F)

set.seed(6233)
fc2plot_sDE_p_fc_bot <- ggplot() +
  geom_point(data = fc2dat %>% filter(padj < 0.05, FC2 > 0.5), aes(KO_Gene, FC2)) +
  geom_text_repel(data = fc2dat %>% filter(padj < 0.05, FC2 > 0.5), aes(KO_Gene, FC2, label = Paralog)) +
  theme_classic() +
  xlab('CRISPR target') +
  ylab('Paralog log2(fold change)') +
  # ggtitle('DESeq2 filters: p.adj < 0.05, LFC > 0.5') +
  geom_hline(yintercept = 0.5, linetype = 2) +
  ylim(c(0, 3))

fc2plot_sDE_p_fc_top <- ggplot() +
  geom_point(data = fc2dat %>% filter(padj < 0.05, FC2 > 0.5), aes(KO_Gene, FC2)) +
  geom_text_repel(data = fc2dat %>% filter(padj < 0.05, FC2 > 0.5), aes(KO_Gene, FC2, label = Paralog)) +
  theme_classic() +
  xlab('CRISPR target') +
  ylab('Paralog log2(fold change)') +
  ggtitle('DESeq2 filters: p.adj < 0.05, LFC > 0.5') +
  scale_y_continuous(breaks = c(8,9,10)) +
  ylim(c(8,10)) 

pdf('de_analysis/summary_plots_drafts/fc2plot_onlySigDE_splitY_p_fc.pdf')
grid.arrange(fc2plot_sDE_p_fc_top, fc2plot_sDE_p_fc_bot, nrow = 2, heights = c(2,3))
dev.off()

svglite('de_analysis/summary_plots_drafts/fc2plot_onlySigDE_splitY_p_fc.svg', width = 11, height = 6)
grid.arrange(fc2plot_sDE_p_fc_top, fc2plot_sDE_p_fc_bot, nrow = 2, heights = c(2,3))
dev.off()

perup_summary_p_fc_nonhit <- perup_summary_p_fc %>% filter(!(KO_Gene %in% fc2dat$KO_Gene))
fc2dat_nh <- list()
corparas_nh <- list()
for (ko in perup_summary_p_fc_nonhit$KO_Gene) {
  
  t_hit <- perup_summary_p_fc_nonhit %>% filter(KO_Gene == ko)
  
  tda <- as_tibble(read.csv(paste0('de_analysis/p_fc/FC2/tabular_format/', t_hit$GEO_ID,'.csv'))) %>%
    filter(KO_Gene == ko) %>%
    mutate(GEO_ID = t_hit$GEO_ID) %>%
    dplyr::select(-X) %>%
    filter(!is.na(padj))
  
  if(is.null(dim(fc2dat_nh))) {
    fc2dat_nh <- tda
  } else {
    fc2dat_nh %<>% bind_rows(tda)
  }
  
  if (length(tda$Perc_ID) > 2) {
    cort <- cor.test(tda$Perc_ID, tda$FC2)
    if(is.null(dim(corparas_nh))){
      corparas_nh <- tibble(KO_Gene = ko,
                         Cor = cort$estimate,
                         P.val = cort$p.value)
    } else {
      corparas_nh %<>% bind_rows(tibble(KO_Gene = ko,
                                     Cor = cort$estimate,
                                     P.val = cort$p.value))
    }
  }
}

write.csv(fc2dat_nh, file ='de_analysis/summary_plots_drafts/foldchanges_nonhitsonly_anypara.csv', quote = F, row.names = F)



para_homology_FC2_cor <- ggplot() + 
  geom_point(data = fc2dat, aes(Perc_ID, FC2)) +
  geom_text(data = corparas, aes(x = 50, y = 5, label = paste0('R = ', as.character(round(Cor,3)), '\nP-val = ', as.character(round(P.val,3))))) +
  facet_wrap(~KO_Gene) +
  theme_classic() +
  xlab('% sequence similarity') +
  ylab('Paralog log2(fold change)')
ggsave(para_homology_FC2_cor, filename = 'de_analysis/summary_plots_drafts/para_homology_FC2_cor.pdf', width = 6, height = 6)
ggsave(para_homology_FC2_cor, filename = 'de_analysis/summary_plots_drafts/para_homology_FC2_cor.svg', width = 6, height = 6)

