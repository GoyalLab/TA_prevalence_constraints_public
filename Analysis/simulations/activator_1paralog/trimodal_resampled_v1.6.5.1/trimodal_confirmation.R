library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
library(grid)
library(diptest)
library(e1071)
library(ggrepel)
library(corrplot)
library(tibble)
library(svglite)
library(entropy)
library(ggalluvial)
library(rgl)
source('~/code/Github/grn_nitc/Functions/grn_analysis_utilities.R')


# edit as needed
datadir51 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.1/samples/'
plotdir51 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.1/exploratory_analysis/'
setwd(datadir51)

if(!dir.exists(plotdir51)){
  dir.create(plotdir51)
}
paramsets51 <- 1:100
lhs_sets51 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets51 %<>%
  mutate(paramset = 1:nrow(lhs_sets51))

lhs_sets51_Hn5<- lhs_sets51 %>% filter(Hill_coefficient_n < 5, paramset <= 9900)

paramsets51 <- lhs_sets51_Hn5$paramset

# calculate stats
allstats51 <- list()
allparams51 <- list()
for (paramset in paramsets51){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 1 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir51, 'distributions_q300_v1.6.5.1_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  if(is.null(dim(allstats51))) {
    allstats51 <- spstats
  } else {
    allstats51 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams51))) {
    allparams51 <- params
  } else {
    allparams51 %<>% bind_rows(params)
  }
  
}

write.csv(tibble(paramset = 1:100), paste0(plotdir51, 'trimodal_sets_manual.csv'), quote = F, row.names = F)

# sample 100 random parameter sets from 1.6.5.2
set.seed(362)
sampsets <- sample(1:10000, 100)
setwd(datadir52)
lhs_sets52 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets52 %<>%
  mutate(paramset = 1:nrow(lhs_sets52))
lhs_sets52_samp <- lhs_sets52 %>%
  filter(paramset %in% sampsets)

for (paramset in sampsets){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 1 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir51, 'distributions_q300_v1.6.5.2_paramset_', as.character(paramset), '.pdf'))
  }
  
}

write.csv(tibble(paramset = sampsets[order(sampsets)]), paste0(plotdir51, 'random_sets_1.6.5.2_manual.csv'), quote = F, row.names = F)

manTrimodalCheck <- bind_rows(
  as_tibble(read.csv(paste0(plotdir51, 'trimodal_sets_manual.csv'))) %>% mutate(version = 'resampled_1.6.5.1'),
  as_tibble(read.csv(paste0(plotdir51, 'random_sets_1.6.5.2_manual.csv'))) %>% mutate(version = 'original_1.6.5.2')
)

trimCheckStackedBar <- ggplot(manTrimodalCheck, aes(version, fill = isTrimodal)) + 
  geom_bar() +
  scale_fill_brewer(palette = 'Dark2') +
  theme_classic()
ggsave(trimCheckStackedBar, file = paste0(plotdir51, 'trimodal_manualcheck_stackedbar.pdf'))

## plot parameter sets on full sample space (2 parameters at a time)

all_sampled_lhs_5152 <- bind_rows(
  lhs_sets51 %>%
    mutate(version = 'v1.6.5.1'),
  lhs_sets52_samp %>%
    mutate(version = 'v1.6.5.2')
)
colors5152 <- c('blue','black')
names(colors5152) <- c('v1.6.5.1', 'v1.6.5.2')

refset <- lhs_sets51[1,]

refset_radii <- c(0.1, 0.01, 0.1, 0.1, 1, 0.1, 0.1, 0.1)
names(refset_radii) <- colnames(lhs_sets51)[1:8]

basalnitc_addonratio_paramscatter <- ggplot(all_sampled_lhs_5152, aes(log10(basal_nitc_on_ratio), log10(A1_Aprime1_addon_ratio), color = version)) +
  geom_point() +
  annotate('rect', 
           xmin = log10(refset$basal_nitc_on_ratio - refset_radii['basal_nitc_on_ratio']), 
           xmax = log10(refset$basal_nitc_on_ratio + refset_radii['basal_nitc_on_ratio']), 
           ymin = log10(refset$A1_Aprime1_addon_ratio - refset_radii['A1_Aprime1_addon_ratio']),
           ymax = log10(refset$A1_Aprime1_addon_ratio + refset_radii['A1_Aprime1_addon_ratio']),
           alpha = 0.1,
           fill = 'blue') +
  theme_classic() +
  scale_color_manual(values = c('v1.6.5.2' = 'black', 'v1.6.5.1' = 'blue')) +
  ggtitle('Trimodal confirmation parameter sets\n100 sets in subspace (blue)\n100 sets in full original space (black)') +
  theme(legend.position = 'none')

prodon_onoff_paramscatter <- ggplot(all_sampled_lhs_5152, aes(log10(A1_Aprime_prodon_ratio), log10(onbasalA1_off_ratio), color = version)) +
  geom_point() +
  annotate('rect', 
           xmin = log10(refset$A1_Aprime_prodon_ratio - refset_radii['A1_Aprime_prodon_ratio']), 
           xmax = log10(refset$A1_Aprime_prodon_ratio + refset_radii['A1_Aprime_prodon_ratio']), 
           ymin = log10(refset$onbasalA1_off_ratio - refset_radii['onbasalA1_off_ratio']),
           ymax = log10(refset$onbasalA1_off_ratio + refset_radii['onbasalA1_off_ratio']),
           alpha = 0.1,
           fill = 'blue') +
  theme_classic() +
  scale_color_manual(values = c('v1.6.5.2' = 'black', 'v1.6.5.1' = 'blue')) +
  ggtitle('Trimodal confirmation parameter sets\n100 sets in subspace (blue)\n100 sets in full original space (black)') +
  theme(legend.position = 'none')

allprodon_A1add_paramscatter <- ggplot(all_sampled_lhs_5152, aes(log10(r_prod_on), log10(r_addon_byA1_B1), color = version)) +
  geom_point() +
  annotate('rect',
           xmin = log10(refset$r_prod_on - refset_radii['r_prod_on']), 
           xmax = log10(refset$r_prod_on + refset_radii['r_prod_on']), 
           ymin = log10(refset$r_addon_byA1_B1 - refset_radii['r_addon_byA1_B1']),
           ymax = log10(refset$r_addon_byA1_B1 + refset_radii['r_addon_byA1_B1']),
           alpha = 0.1,
           fill = 'blue') +
  theme_classic() +
  scale_color_manual(values = c('v1.6.5.2' = 'black', 'v1.6.5.1' = 'blue')) +
  ggtitle('Trimodal confirmation parameter sets\n100 sets in subspace (blue)\n100 sets in full original space (black)') +
  theme(legend.position = 'none')

onbasal_hill_paramscatter <- ggplot(all_sampled_lhs_5152, aes(log10(r_onbasal_A1), log10(Hill_coefficient_n), color = version)) +
  geom_point() +
  annotate('rect',
           xmin = log10(refset$r_onbasal_A1 - refset_radii['r_onbasal_A1']), 
           xmax = log10(refset$r_onbasal_A1 + refset_radii['r_onbasal_A1']), 
           ymin = log10(refset$Hill_coefficient_n - refset_radii['Hill_coefficient_n']),
           ymax = log10(refset$Hill_coefficient_n + refset_radii['Hill_coefficient_n']),
           alpha = 0.1,
           fill = 'blue') +
  theme_classic() +
  scale_color_manual(values = c('v1.6.5.2' = 'black', 'v1.6.5.1' = 'blue')) +
  ggtitle('Trimodal confirmation parameter sets\n100 sets in subspace (blue)\n100 sets in full original space (black)') +
  theme(legend.position = 'none')

ggsave(basalnitc_addonratio_paramscatter, file = paste0(plotdir51, 'subspace_basalnitc_addonratio_paramscatter.pdf'), width = 4, height = 4)
ggsave(prodon_onoff_paramscatter,  file = paste0(plotdir51, 'subspace_prodon_onoff_paramscatter.pdf'), width = 4, height = 4)
ggsave(allprodon_A1add_paramscatter,  file = paste0(plotdir51, 'subspace_allprodon_A1add_paramscatter.pdf'), width = 4, height = 4)
ggsave(onbasal_hill_paramscatter,  file = paste0(plotdir51, 'subspace_onbasal_hill_paramscatter.pdf'), width = 4, height = 4)

# 3d doesn't work flexibly enough   
# basalnitc_addonratio_prodonratio_paramscatter <- plot3d(log10(all_sampled_lhs_5152$basal_nitc_on_ratio),
#                                                         log10(all_sampled_lhs_5152$A1_Aprime1_addon_ratio),
#                                                         log10(all_sampled_lhs_5152$A1_Aprime_prodon_ratio),
#                                                         xlab = 'log10(basal_nitc_on_ratio)',
#                                                         ylab = 'log10(A1_Aprime1_addon_ratio)',
#                                                         zlab = 'log10(A1_Aprime_prodon_ratio)',
#                                                         type = 'p',
#                                                         col = colors5152[all_sampled_lhs_5152$version])
# 
# c3d1 <- cube3d(color = 'blue', alpha = 0.2) %>%
#   translate3d(log10(lhs_sets51$basal_nitc_on_ratio[1]),log10(lhs_sets51$A1_Aprime1_addon_ratio[1]),log10(lhs_sets51$A1_Aprime_prodon_ratio[1])) %>%
#   scale3d(refset_radii['basal_nitc_on_ratio'], refset_radii['A1_Aprime1_addon_ratio'], refset_radii['A1_Aprime_prodon_ratio'])
# shade3d(c3d1)


# pull specific parameter sets and plot histograms and traces for figure (svg files)

setsToPlot <- tibble(
  version = c('1.6.5.1'),
  paramset = c(24),
  classID = c('trimodal')
)

if(!dir.exists(paste0(plotdir51, 'panel_drafts/'))){
  dir.create(paste0(plotdir51, 'panel_drafts/'))
}

for (ind in 1:nrow(setsToPlot)){
  
  pset = setsToPlot$paramset[ind]
  ver = setsToPlot$version[ind]
  classlab = setsToPlot$classID[ind]
  
  if(!dir.exists(paste0(plotdir51, 'panel_drafts/', classlab))){
    dir.create(paste0(plotdir51, 'panel_drafts/', classlab))
  }
  if(!dir.exists(paste0(plotdir51, 'panel_drafts/', classlab, '/if_panel_1C/'))){
    dir.create(paste0(plotdir51, 'panel_drafts/', classlab, '/if_panel_1C/'))
  }  
  if(!dir.exists(paste0(plotdir51, 'panel_drafts/', classlab, '/supp/'))){
    dir.create(paste0(plotdir51, 'panel_drafts/', classlab, '/supp/'))
  }

  tracedir = datadir51
  setwd(tracedir)
  species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0n <- plot_traces_ver_nonons(species, 7500, 8000, 'WT/WT')
  traceplot1n <- plot_traces_ver_nonons(species, 107500, 108000, 'WT/MUT')
  traceplot2n <- plot_traces_ver_nonons(species, 207500, 208000, 'MUT/MUT')
  
  f0n<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_noNonsense.svg')
  f1n<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_noNonsense.svg')
  f2n<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_noNonsense.svg')
  
  ggsave(plot = traceplot0n, f0n, width = 8, height = 6)
  ggsave(plot = traceplot1n, f1n, width = 8, height = 6)
  ggsave(plot = traceplot2n, f2n, width = 8, height = 6)
  
  
  traceplot0 <- plot_traces_ver(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces_ver(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces_ver(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir51, 'panel_drafts/', classlab, '/if_panel_1C/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT.svg')
  f1<-paste0(plotdir51, 'panel_drafts/', classlab, '/if_panel_1C/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT.svg')
  f2<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT.svg')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  
  traceplot0sb <- plot_traces_ver_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1sb <- plot_traces_ver_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2sb <- plot_traces_ver_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0sb<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus_short.svg')
  f1sb<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus_short.svg')
  f2sb<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus_short.svg')
  
  ggsave(plot = traceplot0sb, f0sb, width = 8, height = 6)
  ggsave(plot = traceplot1sb, f1sb, width = 8, height = 6)
  ggsave(plot = traceplot2sb, f2sb, width = 8, height = 6)
  
  traceplot0nb <- plot_traces_ver_nonons_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1nb <- plot_traces_ver_nonons_Bfocus(species, 107800, 108000, 'WT/MUT')
  traceplot2nb <- plot_traces_ver_nonons_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0nb<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_noNonsense_Bfocus_short.svg')
  f1nb<-paste0(plotdir51, 'panel_drafts/', classlab, '/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_noNonsense_Bfocus_short.svg')
  f2nb<-paste0(plotdir51, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_noNonsense_Bfocus_short.svg')
  
  ggsave(plot = traceplot0nb, f0nb, width = 8, height = 6)
  ggsave(plot = traceplot1nb, f1nb, width = 8, height = 6)
  ggsave(plot = traceplot2nb, f2nb, width = 8, height = 6)
  
  
  
  species_sample <- species %>%
    mutate(paramset = pset, version = ver) %>%
    filter( (time+1) %% 300 == 0, (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() +
    geom_rug() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot, file = paste0(plotdir51, 'panel_drafts/', classlab, '/supp/histogram_in_wtmut_distributions_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.svg'))
  
  dist_plot2<-ggplot(species_sample %>% filter(mutated_alleles < 2), aes(abundance)) +
    geom_histogram() +
    geom_rug() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot2, file = paste0(plotdir51, 'panel_drafts/', classlab, '/if_panel_1C/histogram_in_wtmut_distributions_WTWT_WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.svg'))
  
  dist_plot3 <- ggplot(species_sample %>% filter(mutated_alleles == 1, product == 'B1'), aes(abundance)) +
    geom_histogram() +
    geom_rug() +
    ggtitle(paste0('Version ', as.character(ver), ' Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot3, file = paste0(plotdir51, 'panel_drafts/', classlab, '/histogram_in_wtmut_distribution_B1WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.svg'))
  
}

# compare parameter ranges
paramspace_linear51 <- lhs_sets51 %>%
  pivot_longer(!paramset, 
               names_to = 'parameter',
               values_to = 'value') %>%
  group_by(parameter) %>%
  summarise(minimum = min(value),
            maximum = max(value))
paramspace_log51 <- paramspace_linear51 %>%
  mutate(minimum = log10(minimum),
         maximum = log10(maximum),
         range = maximum - minimum)

isbimodal_het_log_paramranges <- as_tibble(read.csv(paste0(plotdir52, 'stats_class_assignment_check_v', as.character(anver),'/isbimodal_het_log_paramranges.csv'), stringsAsFactors = F, header = T))
isunimodalsymmetric_het_log_paramranges <- as_tibble(read.csv(paste0(plotdir52, 'stats_class_assignment_check_v', as.character(anver),'/isunimodalsymmetric_het_log_paramranges.csv'), stringsAsFactors = F, header = T))
isleftskewed_het_log_paramranges <- as_tibble(read.csv(paste0(plotdir52, 'stats_class_assignment_check_v', as.character(anver),'/isleftskewed_het_log_paramranges.csv'), stringsAsFactors = F, header = T))
isrightskewed_het_log_paramranges <- as_tibble(read.csv(paste0(plotdir52, 'stats_class_assignment_check_v', as.character(anver),'/isrightskewed_het_log_paramranges.csv'), stringsAsFactors = F, header = T))
islowaverage_het_log_paramranges <- as_tibble(read.csv(paste0(plotdir52, 'stats_class_assignment_check_v', as.character(anver),'/islowaverage_het_log_paramranges.csv'), stringsAsFactors = F, header = T))

dt_ranges <- bind_rows(isbimodal_het_log_paramranges %>% mutate(shape = 'bimodal'),
                       isunimodalsymmetric_het_log_paramranges %>% mutate(shape = 'unimodal-symmetric')) %>%
  bind_rows(isleftskewed_het_log_paramranges %>% mutate(shape = 'left-skewed')) %>%
  bind_rows(isrightskewed_het_log_paramranges %>% mutate(shape = 'right-skewed')) %>%
  bind_rows(islowaverage_het_log_paramranges %>% mutate(shape = 'low-average'))

max_ranges <- dt_ranges %>%
  group_by(parameter) %>%
  summarise(max_range = max(range))

dt_ranges %<>% left_join(max_ranges, by = 'parameter')

ranges_histogram_includingmax <- ggplot() +
  geom_histogram(data = dt_ranges, aes(range)) +
  geom_vline(data = paramspace_log51, aes(xintercept = range), color = 'blue') +
  facet_wrap(~parameter) +
  theme_classic()

ranges_histogram_excludingmax <- ggplot() +
  geom_histogram(data = dt_ranges %>%
                   filter(range != max_range), aes(range)) +
  geom_vline(data = paramspace_log51, aes(xintercept = range), color = 'blue') +
  facet_wrap(~parameter) +
  theme_classic()

ggsave(ranges_histogram_includingmax, file = paste0(plotdir51, 'panel_drafts/ranges_histogram_includingmax.svg'))
ggsave(ranges_histogram_excludingmax, file = paste0(plotdir51, 'panel_drafts/ranges_histogram_excludingmax.svg'))

# pick minimum plausible subspace (actual ranges from other trees for each parameter, centered around first trimodal set)
larger_samp <- dt_ranges %>%
  group_by(parameter) %>%
  summarise(min_range = min(range)) %>%
  mutate(radius = min_range/2)

orig_set <- lhs_sets51[1,]

larger_resamp_bounds <- as_tibble(orig_set) %>%
  pivot_longer(-paramset, names_to = 'parameter', values_to = 'center') %>%
  dplyr::select(-paramset) %>%
  mutate(log10_center = log10(center)) %>%
  inner_join(larger_samp, by = 'parameter') %>%
  mutate(min_val_log10 = log10_center - radius,
         max_val_log10 = log10_center + radius)
write.csv(larger_resamp_bounds, paste0(plotdir51, 'panel_drafts/larger_resamp_bounds.csv'), quote =F, row.names = F)
larger_resamp_bounds <- as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/larger_resamp_bounds.csv', stringsAsFactors = F))

more_larger_samp <- dt_ranges %>%
  group_by(parameter) %>%
  mutate(max_range = max(range)) %>%
  filter(range != max_range) %>%
  summarise(min_range = min(range),
            perc10_range = quantile(range, 0.1),
            perc50_range = quantile(range, 0.5)) %>%
  mutate(radius_min = min_range/2,
         radius_10 = perc10_range/2,
         radius_50 = perc50_range/2)
more_larger_resamp_bounds <- as_tibble(orig_set) %>%
  pivot_longer(-paramset, names_to = 'parameter', values_to = 'center') %>%
  dplyr::select(-paramset) %>%
  mutate(log10_center = log10(center)) %>%
  inner_join(more_larger_samp, by = 'parameter') %>%
  mutate(min_val_min_log10 = log10_center - radius_min,
         max_val_min_log10 = log10_center + radius_min,
         min_val_10_log10 = log10_center - radius_10,
         max_val_10_log10 = log10_center + radius_10,
         min_val_50_log10 = log10_center - radius_50,
         max_val_50_log10 = log10_center + radius_50)

write.csv(more_larger_resamp_bounds, paste0(plotdir51, 'panel_drafts/more_larger_resamp_bounds.csv'), quote =F, row.names = F)
more_larger_resamp_bounds <- as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/more_larger_resamp_bounds.csv', stringsAsFactors = F))

# edit as needed
datadir512 <- '~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/samples/'
plotdir512 <- '~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/exploratory_analysis/'
setwd(datadir512)

if(!dir.exists(plotdir512)){
  dir.create(plotdir512)
}
paramsets512 <- 1:100
lhs_sets512 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets512 %<>%
  mutate(paramset = 1:nrow(lhs_sets512))

lhs_sets512_Hn5<- lhs_sets512 %>% filter(Hill_coefficient_n < 5, paramset <= 9900)

paramsets512 <- lhs_sets512_Hn5$paramset

# calculate stats
allstats512 <- list()
allparams512 <- list()
for (paramset in paramsets512){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 1 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir512, 'distributions_q300_v1.6.5.1_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  if(is.null(dim(allstats512))) {
    allstats512 <- spstats
  } else {
    allstats512 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams512))) {
    allparams512 <- params
  } else {
    allparams512 %<>% bind_rows(params)
  }
  
}

write.csv(tibble(paramset = 1:100), paste0(plotdir512, 'trimodal_sets_manual.csv'), quote = F, row.names = F)

# edit as needed
datadir513 <- '~/Documents/grn_nitc_data/simulations/v1.6.5.1.3/samples/'
plotdir513 <- '~/Documents/grn_nitc_data/simulations/v1.6.5.1.3/exploratory_analysis/'
setwd(datadir513)

if(!dir.exists(plotdir513)){
  dir.create(plotdir513)
}
paramsets513 <- 1:100
lhs_sets513 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets513 %<>%
  mutate(paramset = 1:nrow(lhs_sets513))

lhs_sets513_Hn5<- lhs_sets513 %>% filter(Hill_coefficient_n < 5, paramset <= 9900)

paramsets513 <- lhs_sets513_Hn5$paramset

# calculate stats
allstats513 <- list()
allparams513 <- list()
for (paramset in paramsets513){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 1 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir513, 'distributions_q300_v1.6.5.1.3_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  if(is.null(dim(allstats513))) {
    allstats513 <- spstats
  } else {
    allstats513 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams513))) {
    allparams513 <- params
  } else {
    allparams513 %<>% bind_rows(params)
  }
  
}

write.csv(tibble(paramset = 1:100), paste0(plotdir513, 'trimodal_sets_manual.csv'), quote = F, row.names = F)


# edit as needed
datadir514 <- '~/Documents/grn_nitc_data/simulations/v1.6.5.1.4/samples/'
plotdir514 <- '~/Documents/grn_nitc_data/simulations/v1.6.5.1.4/exploratory_analysis/'
setwd(datadir514)

if(!dir.exists(plotdir514)){
  dir.create(plotdir514)
}
paramsets514 <- 1:100
lhs_sets514 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets514 %<>%
  mutate(paramset = 1:nrow(lhs_sets514))

lhs_sets514_Hn5<- lhs_sets514 %>% filter(Hill_coefficient_n < 5, paramset <= 9900)

paramsets514 <- lhs_sets514_Hn5$paramset

# calculate stats
allstats514 <- list()
allparams514 <- list()
for (paramset in paramsets514){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 1 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir514, 'distributions_q300_v1.6.5.1.3_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  if(is.null(dim(allstats514))) {
    allstats514 <- spstats
  } else {
    allstats514 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams514))) {
    allparams514 <- params
  } else {
    allparams514 %<>% bind_rows(params)
  }
  
}

write.csv(tibble(paramset = 1:100), paste0(plotdir514, 'trimodal_sets_manual.csv'), quote = F, row.names = F)



manTrimodalCheck2 <- bind_rows(
  as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1/trimodal_sets_manual.csv')) %>% mutate(version = 'resampled_custom_smallest'),
  as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1/random_sets_1.6.5.2_manual.csv')) %>% mutate(version = 'original_1.6.5.2'),
) %>%
  bind_rows(as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/exploratory_analysis/trimodal_sets_manual.csv')) %>% mutate(version = 'resampled_minimum_ranges'))%>%
  bind_rows(as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1.3/exploratory_analysis/trimodal_sets_manual.csv')) %>% mutate(version = 'resampled_perc10_ranges'))%>%
  bind_rows(as_tibble(read.csv('~/Documents/grn_nitc_data/simulations/v1.6.5.1.4/exploratory_analysis/trimodal_sets_manual.csv')) %>% mutate(version = 'resampled_perc50_ranges'))



trimCheckStackedBar_multi <- ggplot(manTrimodalCheck2, aes(version, fill = isTrimodal)) + 
  geom_bar() +
  scale_fill_manual(values = c('TRUE' = 'darkblue', 'FALSE' = 'grey50')) +
  theme_classic()
ggsave(trimCheckStackedBar_multi, file = '~/Documents/grn_nitc_data/simulations/v1.6.5.1/trimodal_manualcheck_stackedbar_multirange.pdf')
ggsave(trimCheckStackedBar_multi, file = '~/Documents/grn_nitc_data/simulations/v1.6.5.1/trimodal_manualcheck_stackedbar_multirange.svg')
