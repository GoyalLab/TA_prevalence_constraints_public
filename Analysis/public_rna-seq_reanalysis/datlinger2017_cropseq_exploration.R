library(tidyverse)
library(ggrepel)
library(magrittr)
library(biomaRt)
# library(GenomicFeatures)

setwd('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/GSE92872/')

plotdir <- '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/'
# load and reformat bulk data
bulkcounts <- as_tibble(read.csv(file = 'GSE92872_CROP-seq_Jurkat_TCR.count_matrix.csv', stringsAsFactors = F))

bulk_coldata <- as.data.frame(bulkcounts[1:3,])
rownames(bulk_coldata) <- bulk_coldata$sample_name
bulk_coldata_tall <- t(bulk_coldata)[2:ncol(bulk_coldata),]

bulk_rowdata <- as.data.frame(bulkcounts[5:nrow(bulkcounts),1])
rownames(bulk_rowdata) <- bulk_rowdata$sample_name
bulk_rowdata %<>% dplyr::rename(gene_name = sample_name)

bulk_countmat <- as.data.frame(bulkcounts[5:nrow(bulkcounts),])
rownames(bulk_countmat) <- bulk_countmat$sample_name
bulk_countmat %<>% dplyr::select(-sample_name) %>% as.matrix()

target_genes <- unique(as.character(as.matrix(bulk_coldata[2,2:ncol(bulk_coldata)]))) #includes CTRL (control)


# load CROP-seq data
sccounts <- as_tibble(read.csv(file = 'GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv', stringsAsFactors = F, header = F))

# pull list of paralogs from ensembl
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl") # version 105, Dec 2021
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
                      values = target_genes,
                      mart = human)

# load "gene lengths", (i.e., union of all Refseq transcript exons/UTRs per gene) from hg38 genes
lengthtbl<- as_tibble(read.csv(file = '~/code/grn_nitc/Resources/hg38.ncbiRefSeq.txLengthPerGene.csv', header = T, stringsAsFactors = F)) %>%
  dplyr::rename(gene_name = gene_id)

# # load T-cell network from Marbach et al. (TO DOWNLOAD) - no longer available unfortunately
# tcellNet = as_tibble(read.table("extractedData/Networks/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/21_heart.txt", sep = "\t", header = F, stringsAsFactors = F)) # change to T-cell
# colnames(tcellNet) = c("TF", "Target", "edgeWt")
# tcellNet = unique(tcellNet$TF)

# Load human regulons from Saez lab 2019 Genome Research paper (BiocManager::install("dorothea"); library dorothea)
# requires R v4.0 or later. I think it may be time
regulons <- dorothea::dorothea_hs

# look for nitc in bulk results
bulk_counts_tall <- as.data.frame(bulk_countmat) %>%
  mutate(gene_name = rownames(bulk_countmat)) %>%
  as_tibble() %>%
  pivot_longer(cols = -gene_name, names_to = 'sampleID', values_to = 'count') 

bulk_counts_tall$count <- as.numeric(as.character(bulk_counts_tall$count))

bulk_totalcounts <- bulk_counts_tall %>%
  group_by(sampleID) %>%
  summarise(total_counts = sum(count))

countplot <- bulk_totalcounts %>% dplyr::select(total_counts) %>% unique() %>% ggplot() + geom_histogram(aes(total_counts)) + theme_classic() + ggtitle('total counts per bulk sample')
ggsave(countplot, file = paste0(plotdir, 'total_counts_per_bulk_sample.pdf'))

bulk_counts_tall %<>% 
  inner_join(bulk_totalcounts, by = 'sampleID') %>%
  mutate(RPM = count*1000000/total_counts) %>%
  inner_join(as.data.frame(bulk_coldata_tall) %>% 
               mutate(sampleID = rownames(bulk_coldata_tall)), by = 'sampleID') %>%
  filter(total_counts > 1e6)

# absolute change in TPM after KO per paralog
bulk_counts_tall %<>%
  inner_join(lengthtbl, by = 'gene_name') %>%
  mutate(RPK = count*1e3/length) %>%
  group_by(gene) %>%
  mutate(sumRPKpm = sum(RPK)/1e6,
         TPM = RPK/sumRPKpm) %>%
  dplyr::select(-c('RPK','sumRPKpm','length'))
         
pseud = 1
bulk_FC_perTarget_rpm <- list()
bulk_FC_perTarget_perParalog_rpm <- list()
bulk_FC_perTarget_tpm <- list()
bulk_FC_perTarget_perParalog_tpm <- list()
for(crispr_target in target_genes) {
  
  paralogs <- geneParaList %>%
    filter(external_gene_name == crispr_target)
  
  bulk_sub <- bulk_counts_tall %>%
    filter(gene %in% c(crispr_target, 'CTRL'),
           gene_name %in% c(crispr_target, paralogs$hsapiens_paralog_associated_gene_name))
  
  bulk_sub_sum_rpm <- bulk_sub %>%
    group_by(gene_name, gene, condition) %>%
    summarise(meanRPM = mean(RPM + pseud),
              sdRPM = sd(RPM + pseud),
              semRPM = sdRPM/sqrt(length(RPM)),
              CRISPR_target = crispr_target) %>%
    group_by(CRISPR_target, gene_name, condition) %>%
    pivot_wider(names_from = gene, values_from = c(meanRPM, sdRPM, semRPM)) %>%
    mutate(meanFC = eval(as.symbol(paste0('meanRPM_', crispr_target)))/meanRPM_CTRL,
           sdFC = meanFC*sqrt((eval(as.symbol(paste0('sdRPM_', crispr_target)))/eval(as.symbol(paste0('meanRPM_', crispr_target))))^2 + 
                                (sdRPM_CTRL/meanRPM_CTRL)^2),
           lfc_RPM = log2(meanFC),
           lfc_RPM_up = log2(meanFC+sdFC),
           lfc_RPM_dn = log2(max((meanFC-sdFC),0.001)))
  
  bulk_summary_FC_rpm <- bulk_sub_sum_rpm %>%
    group_by(condition) %>%
    filter(gene_name != crispr_target) %>%
    summarise(gmeanParalogFC = exp(mean(log(meanFC))),
              log2_gmeanParalogFC = log2(gmeanParalogFC),
              nParalogs = length(meanFC)) %>%
    mutate(CRISPR_target = crispr_target)
  
  bulk_sub_sum_tpm <- bulk_sub %>%
    group_by(gene_name, gene, condition) %>%
    summarise(meanTPM = mean(TPM + pseud),
              sdTPM = sd(TPM + pseud),
              semTPM = sdTPM/sqrt(length(TPM)),
              CRISPR_target = crispr_target) %>%
    group_by(CRISPR_target, gene_name, condition) %>%
    pivot_wider(names_from = gene, values_from = c(meanTPM, sdTPM, semTPM)) %>%
    mutate(meanFC_TPM = eval(as.symbol(paste0('meanTPM_', crispr_target)))/meanTPM_CTRL,
           mean_deltaTPM = eval(as.symbol(paste0('meanTPM_', crispr_target))) - meanTPM_CTRL,
           sdFC_TPM = meanFC_TPM*sqrt((eval(as.symbol(paste0('sdTPM_', crispr_target)))/eval(as.symbol(paste0('meanTPM_', crispr_target))))^2 + 
                                (sdTPM_CTRL/meanTPM_CTRL)^2),
           lfc_TPM = log2(meanFC_TPM),
           lfc_TPM_up = log2(meanFC_TPM+sdFC_TPM),
           lfc_TPM_dn = log2(max((meanFC_TPM-sdFC_TPM),0.001)))
  
  bulk_summary_FC_tpm <- bulk_sub_sum_tpm %>%
    group_by(condition) %>%
    filter(gene_name != crispr_target) %>%
    summarise(gmeanParalogFC_TPM = exp(mean(log(meanFC_TPM))),
              mean_mean_deltaTPM = mean(mean_deltaTPM),
              log2_gmeanParalogFC_TPM = log2(gmeanParalogFC_TPM),
              nParalogs = length(meanFC_TPM)) %>%
    mutate(CRISPR_target = crispr_target)
  
  
  p1 <- ggplot(bulk_sub, aes(gene, RPM)) +
    facet_grid(condition~gene_name) +
    geom_point(aes(color = gene)) +
    theme_bw() +
    ggtitle(paste0('bulk gene expression after ', crispr_target, ' knockout'))
  
  ggsave(p1, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO.pdf'), height = 4, width = 2+nrow(paralogs))
  p2 <- ggplot(bulk_sub, aes(gene, log2(RPM+1))) +
    facet_grid(condition~gene_name) +
    geom_point(aes(color = gene)) +
    theme_bw() +
    ylab('log2(RPM+1)') +
    ggtitle(paste0('bulk gene expression after ', crispr_target, ' knockout'))
  
  ggsave(p2, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log.pdf'), height = 4, width = 2+nrow(paralogs))
  
  p3 <- ggplot() +
    geom_point(data = bulk_sub_sum_rpm, aes(log2(meanRPM_CTRL), lfc_RPM)) +
    geom_errorbar(data = bulk_sub_sum_rpm, aes(log2(meanRPM_CTRL), ymin = lfc_RPM_dn, ymax = lfc_RPM_up)) +
    theme_classic() +
    geom_text_repel(data = bulk_sub_sum_rpm, aes(log2(meanRPM_CTRL), lfc_RPM, label = gene_name)) +
    facet_wrap(~condition) +
    ggtitle(paste0('bulk gene expression change relative to control\nafter ', crispr_target, ' knockout')) +
    ylab('log2(fold change) relative to control') +
    xlab('log2(meanRPM+1) in controls')
  
  ggsave(p3, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log_FC.pdf'), height = 4, width = 8)
  
  p4 <- ggplot() +
    geom_point(data = bulk_sub_sum_rpm, aes(log2(meanRPM_CTRL), lfc_RPM)) +
    # geom_errorbar(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), ymin = lfc_RPM_dn, ymax = lfc_RPM_up)) +
    theme_classic() +
    geom_text_repel(data = bulk_sub_sum_rpm, aes(log2(meanRPM_CTRL), lfc_RPM, label = gene_name)) +
    facet_wrap(~condition) +
    ggtitle(paste0('bulk gene expression change relative to control\nafter ', crispr_target, ' knockout')) +
    ylab('log2(fold change) relative to control') +
    xlab('log2(meanRPM+1) in controls')
  
  ggsave(p4, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log_FCnoEB.pdf'), height = 4, width = 8)
  
  p5 <- ggplot() + 
    geom_point(data = bulk_sub_sum_tpm, aes(log2(meanTPM_CTRL), mean_deltaTPM)) +
    geom_text_repel(data = bulk_sub_sum_tpm, aes(log2(meanTPM_CTRL), mean_deltaTPM, label = gene_name), seed = 362) +
    facet_grid(condition ~ CRISPR_target) +
    theme_classic() +
    ggtitle('Absolute increase in paralog expression\nvs paralog average in CTRL')
  ggsave(p5, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log_deltaTPMvsParalogMean.pdf'), height = 4, width = 8)
  
  # p6 <- ggplot() + 
  #   geom_point(data = bulk_sub_sum_tpm, aes(CRISPR_target, mean_deltaTPM)) +
  #   geom_text_repel(data = bulk_sub_sum_tpm, aes(log2(meanTPM_CTRL), mean_deltaTPM, label = gene_name), seed = 362) +
  #   facet_grid(condition ~ CRISPR_target) +
  #   theme_classic() +
  #   ggtitle('Absolute increase in paralog expression\nvs paralog average in CTRL')
  # ggsave(p6, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log_deltaTPMvsTargetMean.pdf'), height = 4, width = 8)
  
  if(is.null(dim(bulk_FC_perTarget_rpm))) {
    bulk_FC_perTarget_rpm <- bulk_summary_FC_rpm
  } else {
    bulk_FC_perTarget_rpm %<>% bind_rows(bulk_summary_FC_rpm)
  }
  
  if(is.null(dim(bulk_FC_perTarget_perParalog_rpm))) {
    bulk_FC_perTarget_perParalog_rpm <- bulk_sub_sum_rpm
  } else {
    bulk_FC_perTarget_perParalog_rpm %<>% bind_rows(bulk_sub_sum_rpm)
  }
  
  if(is.null(dim(bulk_FC_perTarget_tpm))) {
    bulk_FC_perTarget_tpm <- bulk_summary_FC_tpm
  } else {
    bulk_FC_perTarget_tpm %<>% bind_rows(bulk_summary_FC_tpm)
  }
  
  if(is.null(dim(bulk_FC_perTarget_perParalog_tpm))) {
    bulk_FC_perTarget_perParalog_tpm <- bulk_sub_sum_tpm
  } else {
    bulk_FC_perTarget_perParalog_tpm %<>% bind_rows(bulk_sub_sum_tpm)
  }
  
}

bulk_FC_EGR1_stim_FCplot <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog_rpmerParalog %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1', condition == 'stimulated'), aes(CRISPR_target, lfc_RPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_text_repel(data = bulk_FC_perTarget_perParalog_rpmaralog %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1', condition == 'stimulated'), aes(CRISPR_target, lfc_RPM, label = gene_name), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_bar(data = bulk_FC_perTarget_rpm %>% filter(CRISPR_target == 'EGR1', condition == 'stimulated'), aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Log2(fold-change) after KO in points\nLog2(Geometric mean fold change) in bar') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_FC_EGR1_stim_FCplot, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_EGR1_stim_KO_FCplot.pdf', height = 4, width = 2)

bulk_FC_EGR1_FCplot <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog_rpm %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1'), aes(CRISPR_target, lfc_RPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_text_repel(data = bulk_FC_perTarget_perParalog_rpm %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1'), aes(CRISPR_target, lfc_RPM, label = gene_name), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_bar(data = bulk_FC_perTarget_rpmrTarget %>% filter(CRISPR_target == 'EGR1'), aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Log2(fold-change) after KO in points\nLog2(Geometric mean fold change) in bar') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_FC_EGR1_FCplot, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_EGR1_KO_FCplot.pdf', height = 7, width = 2)

bulk_FC_perTarget_FCplot <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog_rpm %>% filter(CRISPR_target != gene_name), aes(CRISPR_target, lfc_RPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_bar(data = bulk_FC_perTarget_rpm, aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Log2(fold-change) after KO in points\nLog2(Geometric mean fold change) in bar') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_FC_perTarget_FCplot, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_KO_log_FCnoEB.pdf', height = 7, width = 7)

bulk_deltaTPM_perTarget <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog_tpm %>% filter(CRISPR_target != gene_name), aes(CRISPR_target, mean_deltaTPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  # geom_bar(data = bulk_FC_perTarget, aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Absolute change in TPM after KO per paralog\nfor each CRISPR target') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (KO - CTRL) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_deltaTPM_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_KO_deltaTPM.pdf', height = 7, width = 7)

bulk_logdeltaTPM_perTarget <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog_tpm %>% filter(CRISPR_target != gene_name), aes(CRISPR_target, log2(mean_deltaTPM)), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  # geom_bar(data = bulk_FC_perTarget, aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Absolute change in TPM after KO per paralog\nfor each CRISPR target') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression log2(KO - CTRL) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_logdeltaTPM_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_KO_logdeltaTPM.pdf', height = 7, width = 7)


paralogs %>% as_tibble()


bulk_FCvsPercID_EGR1_stim_perTarget <- ggplot(bulk_FC_perTarget_perParalog %>% 
                                                filter(CRISPR_target != gene_name,
                                                       CRISPR_target == 'EGR1',
                                                       condition == 'stimulated') %>%
                                                inner_join(geneParaList %>%
                                                             dplyr::rename(CRISPR_target = external_gene_name,
                                                                           gene_name = hsapiens_paralog_associated_gene_name), by = c('CRISPR_target', 'gene_name')), 
                                              aes(hsapiens_paralog_perc_id, lfc_RPM)) +
  geom_point() +
  geom_text_repel(aes(label = gene_name)) +
  facet_wrap(~CRISPR_target) +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab('Log2(fold-change) after KO') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nVs. Percent sequence similarity to reference gene')
ggsave(bulk_FCvsPercID_EGR1_stim_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulk_FCvsPercID_EGR1_stim_perTarget.pdf', height = 4, width = 4)

bulk_FCvsPercID_stim_perTarget <- ggplot() +
  geom_point(data = bulk_FC_perTarget_perParalog %>% 
               filter(CRISPR_target != gene_name,
                      condition == 'stimulated') %>%
               inner_join(geneParaList %>%
                            dplyr::rename(CRISPR_target = external_gene_name,
                                          gene_name = hsapiens_paralog_associated_gene_name), by = c('CRISPR_target', 'gene_name')), 
             aes(hsapiens_paralog_perc_id, lfc_RPM)) +
  facet_wrap(~CRISPR_target) +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab('Log2(fold-change) after KO') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nVs. Percent sequence similarity to reference gene\nStimulated condition')
ggsave(bulk_FCvsPercID_stim_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulk_FCvsPercID_stim_perTarget.pdf', height = 7, width = 7)

bulk_FCvsPercID_unstim_perTarget <- ggplot() +
  geom_point(data = bulk_FC_perTarget_perParalog %>% 
               filter(CRISPR_target != gene_name,
                      condition == 'unstimulated') %>%
               inner_join(geneParaList %>%
                            dplyr::rename(CRISPR_target = external_gene_name,
                                          gene_name = hsapiens_paralog_associated_gene_name), by = c('CRISPR_target', 'gene_name')), 
             aes(hsapiens_paralog_perc_id, lfc_RPM)) +
  facet_wrap(~CRISPR_target) +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab('Log2(fold-change) after KO') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nVs. Percent sequence similarity to reference gene\nUnstimulated condition')
ggsave(bulk_FCvsPercID_unstim_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulk_FCvsPercID_unstim_perTarget.pdf', height = 7, width = 7)

# filter to paralogs that are overexpressed after reference gene KO in bulk
# arbitrarily pick LFC>0.5 for RPM data

upreg_paralogs_bulk_stim <- bulk_FC_perTarget_perParalog_rpm %>%
  filter(lfc_RPM > 0.5, meanRPM_CTRL > 5, condition == 'stimulated') %>%
  dplyr::select(gene_name, CRISPR_target, condition, lfc_RPM, meanRPM_CTRL)

write.csv(upreg_paralogs_bulk_stim, file = '~/code/grn_nitc/datlinger2017_cropseq/upreg_paralogs_bulk_stim.csv', quote = F, row.names = F)

targselfids <- list()
for (targ in target_genes) {
  targselfid <- geneParaList %>%
    dplyr::rename(gene_name = hsapiens_paralog_associated_gene_name,
                  CRISPR_target = external_gene_name) %>%
    filter(CRISPR_target == targ) %>%
    dplyr::select(CRISPR_target, ensembl_gene_id) %>% 
    unique() %>%
    mutate(gene_name = CRISPR_target, hsapiens_paralog_associated_gene_name = ensembl_gene_id)
  
  if(is.null(dim(targselfids))) {
    targselfids <- targselfid
  } else {
    targselfids %<>% bind_rows(targselfid)
  }
  
}

paralogTFcheck <- geneParaList %>%
  dplyr::rename(gene_name = hsapiens_paralog_associated_gene_name,
                CRISPR_target = external_gene_name) %>%
  bind_rows(targselfids) %>%
  mutate(CRISPR_target_isTF = ensembl_gene_id %in% tfTab$ensembl_gene_id,
         paralog_isTF = hsapiens_paralog_ensembl_gene %in% tfTab$ensembl_gene_id
  ) %>% as_tibble() %>% 
  dplyr::select(CRISPR_target, gene_name, CRISPR_target_isTF, paralog_isTF)

# draft map of reference genes and paralogs to their downstream targets, if TFs, and if in Marbach T-cell network
paralog_and_target_regulons <- paralogTFcheck %>%
  dplyr::rename(tf = gene_name) %>%
  left_join(regulons %>% filter(confidence %in% c('A', 'B', 'C')), by = 'tf')

# process single-cell data
# sccounts$condition[1:10]

sccounts_coldata <- sccounts[1:5,] %>% as.data.frame()
rownames(sccounts_coldata) <- sccounts_coldata[,1]
sccounts_coldata_tall <- t(sccounts_coldata)
sccounts_coldata_tall <- sccounts_coldata_tall[2:nrow(sccounts_coldata_tall),]

sccounts <- sccounts[7:nrow(sccounts),]
sccounts %<>% as.data.frame()
rownames(sccounts) <- sccounts$V1
rodat_sccounts <- sccounts$V1
colnames(sccounts) <- c('GENE', sccounts_coldata_tall[,'cell'])

sccounts <- sccounts[,2:ncol(sccounts)]
sccounts <- sapply(sccounts, as.numeric)

sctotals <- colSums(sccounts) # total counts per cell

scrpms <- apply(sccounts, 2 , function(x){x*1e6/sum(x)}) # rpms

rownames(scrpms) <- rodat_sccounts
# scrpms %<>% as_tibble() %>% mutate(gene_name = rodat_sccounts)

minCountsPerCell = 1000
# next: separate by stimulated vs unstimulated; start w stimulated 
for (crispr_target in target_genes) {
  
  control_cells <- sccounts_coldata_tall %>%
    as_tibble() %>%
    filter(gene == 'CTRL',
           condition == 'stimulated')
  
  targeted_cells <- sccounts_coldata_tall %>%
    as_tibble() %>%
    filter(gene == crispr_target,
           condition == 'stimulated')
  
  paralogs_of_target <- geneParaList %>%
    filter(external_gene_name == crispr_target,
           hsapiens_paralog_associated_gene_name %in% rownames(scrpms))
  
  if(nrow(paralogs_of_target) > 0) {
    
    upreg_paralogs_bulk_stim_targ <- paralogs_of_target %>%
      dplyr::rename(gene_name = hsapiens_paralog_associated_gene_name,
                    CRISPR_target = external_gene_name) %>%
      left_join(upreg_paralogs_bulk_stim, by = c('gene_name', 'CRISPR_target'))
    
    control_rpms <- scrpms[c(paralogs_of_target$hsapiens_paralog_associated_gene_name, crispr_target), control_cells$cell]
    targeted_rpms <- scrpms[c(paralogs_of_target$hsapiens_paralog_associated_gene_name, crispr_target), targeted_cells$cell]
    
    control_rpms %<>% as.data.frame()
    control_rpms$gene_name = rownames(control_rpms)
    
    targeted_rpms %<>% as.data.frame()
    targeted_rpms$gene_name = rownames(targeted_rpms)
    
    control_rpms_tall <- control_rpms %>%
      pivot_longer(cols = -gene_name, names_to = 'cell', values_to = 'rpm') %>%
      mutate(CRISPR_target = 'CTRL')
    targeted_rpms_tall <- targeted_rpms %>%
      pivot_longer(cols = -gene_name, names_to = 'cell', values_to = 'rpm') %>%
      mutate(CRISPR_target = crispr_target)
    
    max_exp <- max(c(control_rpms_tall$rpm, targeted_rpms_tall$rpm))
    
    sc_histograms <- ggplot() +
      geom_histogram(data = bind_rows(control_rpms_tall, targeted_rpms_tall), aes(rpm)) +
      geom_text(data = upreg_paralogs_bulk_stim_targ, aes(x=max_exp/2, y = 100, label = ifelse(!is.na(meanRPM_CTRL),
                                                                                               paste0('Bulk mean ', as.character(meanRPM_CTRL),'\nLFC = ', as.character(lfc_RPM)),''))) +
      facet_grid(CRISPR_target ~ gene_name, scales = 'free') +
      theme_classic() +
      ggtitle(paste0('Paralog expression before and after ', crispr_target ,' mutation\n', as.character(nrow(control_cells)), ' control cells, ', as.character(nrow(targeted_cells)), ' targeted cells'))
    ggsave(sc_histograms, file = paste0(plotdir, 'singleCell_', crispr_target, '_paralog_stimulated_histograms.pdf'), width = 2+0.8*nrow(paralogs_of_target), height = 7)
    
    sc_densityplots <- ggplot() +
      geom_histogram(data = bind_rows(control_rpms_tall, targeted_rpms_tall), aes(x=rpm, y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
      geom_text(data = upreg_paralogs_bulk_stim_targ, aes(x=max_exp/2, y = 0.8, label = ifelse(!is.na(meanRPM_CTRL),
                                                                                               paste0('Bulk mean ', as.character(meanRPM_CTRL),'\nLFC = ', as.character(lfc_RPM)),''))) +
      facet_grid(CRISPR_target ~ gene_name, scales = 'free') +
      theme_classic() +
      ggtitle(paste0('Paralog expression before and after ', crispr_target ,' mutation\n', as.character(nrow(control_cells)), ' control cells, ', as.character(nrow(targeted_cells)), ' targeted cells'))
    ggsave(sc_densityplots, file = paste0(plotdir, 'singleCell_', crispr_target, '_paralog_stimulated_densityplots.pdf'), width = 2+0.8*nrow(paralogs_of_target), height = 7)
    
    sc_violins <- ggplot() +
      geom_violin(data = bind_rows(control_rpms_tall, targeted_rpms_tall), aes(x=gene_name, y=rpm, fill=CRISPR_target, color=CRISPR_target)) +
      geom_text(data = upreg_paralogs_bulk_stim_targ, aes(x=gene_name, y = max_exp*0.9, label = ifelse(!is.na(meanRPM_CTRL),
                                                                                               paste0('Bulk mean ', as.character(meanRPM_CTRL),'\nLFC = ', as.character(lfc_RPM)),''))) +
      # facet_grid(CRISPR_target ~ gene_name, scales = 'free') +
      theme_classic() +
      ggtitle(paste0('Paralog expression before and after ', crispr_target ,' mutation\n', as.character(nrow(control_cells)), ' control cells, ', as.character(nrow(targeted_cells)), ' targeted cells'))
    ggsave(sc_violins, file = paste0(plotdir, 'singleCell_', crispr_target, '_paralog_stimulated_rpm_violins.pdf'), width = 2+0.8*nrow(paralogs_of_target), height = 7)
    
    sc_violins_log <- ggplot() +
      geom_violin(data = bind_rows(control_rpms_tall, targeted_rpms_tall), aes(x=gene_name, y=log2(rpm+1), fill=CRISPR_target, color=CRISPR_target)) +
      geom_text(data = upreg_paralogs_bulk_stim_targ, aes(x=gene_name, y = log2(max_exp+1)*0.9, label = ifelse(!is.na(meanRPM_CTRL),
                                                                                                       paste0('Bulk mean ', as.character(meanRPM_CTRL),'\nLFC = ', as.character(lfc_RPM)),''))) +
      # facet_grid(CRISPR_target ~ gene_name, scales = 'free') +
      theme_classic() +
      ggtitle(paste0('Paralog expression (log2(RPM+1)) before and after ', crispr_target ,' mutation\n', as.character(nrow(control_cells)), ' control cells, ', as.character(nrow(targeted_cells)), ' targeted cells'))
    ggsave(sc_violins_log, file = paste0(plotdir, 'singleCell_', crispr_target, '_paralog_stimulated_logrpm_violins.pdf'), width = 2+0.8*nrow(paralogs_of_target), height = 7)
    
  }
}


