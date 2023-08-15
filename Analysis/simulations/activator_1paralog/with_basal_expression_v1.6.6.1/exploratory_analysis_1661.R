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
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')


# edit as needed
datadir61 <- '/Volumes/IAMYG2/grn_nitc_data/v1.6.6.1/samples/'
plotdir61 <- '/Volumes/IAMYG2/grn_nitc_data/v1.6.6.1/exploratory_analysis/'
setwd(datadir61)

if(!dir.exists(plotdir61)){
  dir.create(plotdir61)
}
paramsets61 <- 1:9900
lhs_sets61 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets61 %<>%
  mutate(paramset = 1:nrow(lhs_sets61))

# calculate stats
allstats61 <- list()
allparams61 <- list()
for (paramset in paramsets61){
  
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
  
  if(paramset %% 100 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir61, 'distributions_q300_v1.6.6.1_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              skewness = skewness(abundance + 0.01),
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
  
  if(is.null(dim(allstats61))) {
    allstats61 <- spstats
  } else {
    allstats61 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams61))) {
    allparams61 <- params
  } else {
    allparams61 %<>% bind_rows(params)
  }
  
}
write.csv(allstats61 %>% mutate(version = '1.6.6.1'), file = paste0(plotdir61, 'summary_stats.csv'), row.names = F)


# temp: collate all data
setwd(datadir61)
all_species_q300_61 <- list()
for (paramset in paramsets61){
  
  if(paramset %% 100 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
  }
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset,
           version = '1.6.6.1') %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, version, time, mutated_alleles)
  
  
  if(is.null(dim(all_species_q300_61))) {
    all_species_q300_61 <- species_sample
  } else {
    all_species_q300_61 %<>% bind_rows(species_sample)
  }
  
}

write.csv(all_species_q300_61, file = paste0('/Volumes/IAMYG2/grn_nitc_data/v1.6.6.1/all_species_q300.csv'))


# summary stats draft1
allstats_full61 <- allstats61 %>% mutate(version = '1.6.6.1')
allparams_full61 <- allparams61 %>% mutate(paramset = 1:9900, version = '1.6.6.1')
lhs_sets_full61 <- lhs_sets61 %>% mutate(version = '1.6.6.1')

pseud = 0.01

allstats_full61 %<>% mutate(skewness = ifelse(is.na(skewness), 0, skewness))

compared_stats61 <- allstats_full61 %>% 
  group_by(version, paramset, product, mutated_alleles) %>% 
  pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:entropy90) %>% 
  pivot_wider(names_from = mutated_alleles, values_from = value) %>% 
  mutate(lfc10 = log2((`1`+pseud)/(`0` + pseud)), 
         delta10 = `1`-`0`,
         lfc21 = log2((`2`+pseud)/(`1` + pseud)), 
         delta21 = `2`-`1`,
         lfc20 = log2((`2`+pseud)/(`0` + pseud)), 
         delta20 = `2`-`0`) %>%
  dplyr::select(-c(`0`:`2`)) %>% 
  pivot_longer(names_to = 'compare', values_to = 'diff', cols = lfc10:delta20) %>%
  inner_join(allstats_full61 %>% 
               dplyr::select(mutated_alleles, product, paramset, version, mean_product) %>%
               pivot_wider(names_from = mutated_alleles, values_from = mean_product), by = c('product', 'paramset', 'version')) %>%
  mutate(mean_denom = case_when(
    compare %in% c('lfc10', 'delta10', 'lfc20', 'delta20') ~ `0`,
    compare %in% c('lfc21', 'delta21') ~ `1`))

# filter to Hill n < 5 - now already done for allstats_full1
#compared_stats %<>% ungroup() %>% inner_join(lhs_sets_fall %>% filter(Hill_coefficient_n < 5) %>% dplyr::select(version, paramset), by = c('version','paramset'))

# plot wt/mut summary stats against mean expression

unistats61<-unique(compared_stats61$stat)

for (st in unistats61) {
  
  pvs1 <- ggplot(allstats_full61 %>% inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n<5)) + 
    geom_point(aes(mean_product, eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = 10), linetype = 2) +
    # geom_text(aes(mean_product, eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic() +
    ggtitle(paste0(st, ' vs mean'))
  
  pvs2 <- ggplot(allstats_full61 %>% inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n<5)) + 
    geom_point(aes(log(mean_product), eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = log(10)), linetype = 2) +
    # geom_text(aes(log(mean_product), eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic()  +
    ggtitle(paste0(st, ' vs mean'))
  
  pvs3 <- ggplot(allstats_full61 %>% 
                   inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>% 
                   filter(product == 'B1', Hill_coefficient_n<5, mutated_alleles == 1) %>% 
                   mutate(is10 = mean_product > 10)) + 
    geom_point(aes(log(mean_product), eval(as.symbol(st)), color = is10), alpha = 0.3, stroke = 0) +
    geom_vline(aes(xintercept = log(10)), linetype = 2) +
    scale_color_manual(values = c('grey50', 'black')) +
    # geom_text(aes(log(mean_product), eval(as.symbol(st)), label = as.character(paramset))) +
    # facet_grid(~mutated_alleles) +
    theme_classic()  +
    theme(legend.position = 'none') +
    ylab(st) +
    ggtitle(paste0(st, ' vs mean'))
  
  ggsave(pvs1, file = paste0(plotdir61, 'PerGenotype_', st, '_vs_mean.pdf'), width = 16, height = 8) 
  ggsave(pvs2, file = paste0(plotdir61, 'PerGenotype_', st, '_vs_logmean.pdf'), width = 16, height = 8) 
  ggsave(pvs3, file = paste0(plotdir61, 'WTMUT_', st, '_vs_logmean.pdf'), width = 4, height = 4) 
}


# LOESS
loess_fitted_allstats_all61 <- allstats_full61
for (stat in unistats61[unistats61 != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  statdat <- list()
  statdat1 <- list()
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    # for (ma in 0:2) {
    
    tempdat <- allstats_full61 %>% 
      filter(#mutated_alleles == ma,
        product == gene) %>%
      dplyr::select(mutated_alleles, product, version, paramset, mean_product, stat)
    
    loess1 <- loess(eval(as.symbol(stat)) ~ mean_product, data = tempdat, span = 0.1)
    
    l1dat <- data.frame(mean_product = loess1$x,
                        stat = loess1$fitted,
                        resid = loess1$residuals,
                        version = tempdat$version,
                        paramset = tempdat$paramset,
                        product = gene,
                        mutated_alleles = tempdat$mutated_alleles)
    colnames(l1dat)[2] <- paste0(stat,'_fitted')
    colnames(l1dat)[3] <- paste0(stat,'_residual')
    
    lplot1 <- ggplot() +
      geom_point(data = tempdat, aes(mean_product, eval(as.symbol(stat))), alpha = 0.1) +
      geom_point(data = l1dat, aes(mean_product, eval(as.symbol(paste0(stat,'_fitted')))), color = 'red') +
      theme_classic() +
      ylab(stat) +
      xlab('Mean') +
      ggtitle(paste0(stat, ' vs mean, with LOESS fit to mean\nGene product: ', gene))#, ', mutated alleles: ', as.character(ma)))
    
    lplot2 <- ggplot() +
      geom_point(data = tempdat, aes(log(mean_product), eval(as.symbol(stat))), alpha = 0.1) +
      geom_point(data = l1dat, aes(log(mean_product), eval(as.symbol(paste0(stat,'_fitted')))), color = 'red') +
      theme_classic() +
      ylab(stat) +
      xlab('Log(Mean)') +
      ggtitle(paste0(stat, ' vs log(mean), with LOESS fit to mean\nGene product: ', gene))#, ', mutated alleles: ', as.character(ma)))
    
    # ggsave(lplot1, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.6.1.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_v1.6.2only.pdf'), width = 5, height = 5)
    # ggsave(lplot2, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.6.1.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_log_v1.6.2only.pdf'), width = 5, height = 5)
    
    l1dat$version <- as.character(l1dat$version)
    l1dat$product <- as.character(l1dat$product)
    l1dat$paramset <- as.numeric(l1dat$paramset)
    
    cat('sliding window normalizing...\n') # do this per-gene over all genotypes, rather than all genes over all genotypes...
    l1dat1 <- sliding_window_normalize(as_tibble(l1dat) %>% filter(mean_product>10), 'mean_product', paste0(stat,'_residual'), 50)
    
    if(is.null(dim(statdat))){
      statdat <- l1dat
      statdat1 <- l1dat1
    } else {
      statdat %<>% bind_rows(l1dat)
      statdat1 %<>% bind_rows(l1dat1)
    }
    
    
  }
  
  # Need to split analysis over all genes #
  
  
  loess_fitted_allstats_all61 %<>% left_join(as_tibble(statdat) %>% dplyr::select(-mean_product), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  # 
  # cat('sliding window normalizing...\n') # do this per-gene over all genotypes, rather than all genes over all genotypes...
  # statdat1 <- sliding_window_normalize(as_tibble(statdat) %>% filter(mean_product>10), 'mean_product', paste0(stat,'_residual'), 50)
  # 
  
  loess_fitted_allstats_all61 %<>% left_join(statdat1 %>% dplyr::select(-c('mean_product', paste0(stat,'_residual'), paste0(stat,'_fitted'))), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  
  cat(paste0('Done with ', stat, '\n'))
  
}

write.csv(loess_fitted_allstats_all61, file = paste0(plotdir61, 'loess_fitted_allstats_all_snwRadius100.csv'), quote = F, row.names = F)

# load as needed
# loess_fitted_allstats_all61 <- as_tibble(read.csv(paste0(plotdir61, 'loess_fitted_allstats_all_snwRadius100.csv'), header=T, stringsAsFactors = F))

for (stat in unistats61[unistats61 != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  # statdat <- list()
  
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    
    statdatA = loess_fitted_allstats_all61 %>%
      dplyr::select(mutated_alleles, product, version, paramset, mean_product, stat, as.symbol(paste0(stat, '_residual')), as.symbol(paste0(stat, '_residual_swn'))) %>%
      filter(product == gene, mean_product>10) 
    
    lplot_all_stat_con <- ggplot() +
      geom_point(data = statdatA, aes(log(mean_product), eval(as.symbol(stat))), stroke=0, alpha = 0.05) +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 0, mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'blue') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 1, mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'red') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 2, mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'green') +
      theme_classic()
    
    lplot_all_statLOESS_con <-  ggplot() +
      geom_point(data = statdatA, aes(log(mean_product), eval(as.symbol(paste0(stat,'_residual')))), stroke=0, alpha = 0.05) +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 0,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'blue') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 1,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'red') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 2,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'green') +
      theme_classic()
    
    lplot_all_statLOESSSWN_con <-  ggplot() +
      geom_point(data = statdatA, aes(log(mean_product), eval(as.symbol(paste0(stat,'_residual_swn')))), stroke=0, alpha = 0.05) +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 0,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'blue') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 1,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'red') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 2,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'green') +
      theme_classic()
    
    pdf(paste0(paste0(plotdir61, 'LOESSplots_', stat, 'vsLogMean_', gene,'_v1.6.6.1.pdf')), width = 10, height = 7)
    grid.arrange(lplot_all_stat_con,lplot_all_statLOESS_con,lplot_all_statLOESSSWN_con, ncol=3,
                 top = textGrob(paste0(stat, ' vs. log(mean_product), ', gene, '\nStat, LOESS residual, Squeezed LOESS residual (radius=50)'),gp=gpar(fontsize=20,font=3)))
    dev.off()
    
  }
}

# classification
# focus on LOESS residuals except when specifically indicated (e.g., skewness for exponential dist assignment). sliding window is not normalizing stably enough as intended.

anver <- 5 # increase minimum bimodality_residual filter and make left-skew filter more stringent

bimfilt <- 0.1

entfilt <- 0.15

basic_class_assignment_all61 <- loess_fitted_allstats_all61 %>%
  mutate(class_assignment = case_when(
    mean_product < 10 ~ 'low-average',
    bimodality_coef_residual > bimfilt & bimodality_coef > 0.555 ~ 'bimodal',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & abs(skewness) < 1 ~ 'unimodal symmetric',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness >= 1 ~ 'right-skewed unimodal',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness <= -1 ~ 'left-skewed unimodal'
    
  )) 

basic_class_assignment_all_forSankey61 <- basic_class_assignment_all61 %>%
  dplyr::select(version, paramset, product, mutated_alleles, class_assignment) %>%
  group_by(version, paramset, product) %>%
  pivot_wider(names_from = mutated_alleles, values_from = class_assignment) %>%
  group_by(product, `0`, `1`, `2`) %>%
  summarise(Freq = length(`0`)) %>%
  ungroup() %>%
  mutate(alluvID = 1:length(product)) %>%
  pivot_longer(`0`:`2`, names_to = 'mutated_alleles', values_to = 'class_assignment')

# sample paramsets from the sankey flow to visually inspect accuracy of assignments/changes
set.seed(73245)
basic_class_assignment_all_forSankey_forsamples61 <- basic_class_assignment_all61 %>%
  dplyr::select(version, paramset, product, mutated_alleles, class_assignment)  

classes61 = unique(basic_class_assignment_all61$class_assignment)

if(!dir.exists(paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver)))) {
  dir.create(paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver)))
}

classes_sankey <- ggplot(basic_class_assignment_all_forSankey61, aes(x = mutated_alleles, y=Freq,
                                                                      stratum = class_assignment, alluvium = alluvID, fill = class_assignment, label = class_assignment)) +
  facet_grid(product~.) +
  geom_flow() +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = 'stratum', size = 3) + 
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle('Class assignments before and after mutations\nPositive regulation, log-sampled parameters\nWith basal paralog expression') +
  xlab('Mutated alleles') +
  ylab('Number of parameter sets')
ggsave(classes_sankey, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/classes_sankey.pdf'))

classes_sankey61_B1forfig <- ggplot(basic_class_assignment_all_forSankey61 %>%
                                      filter(product == 'B1',
                                             mutated_alleles %in% c(0,1)), aes(x = mutated_alleles, y=Freq,
                                                                               stratum = class_assignment, alluvium = alluvID, fill = class_assignment, label = class_assignment)) +
  facet_grid(product~.) +
  geom_flow() +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = 'stratum', size = 3) + 
  scale_fill_brewer(palette = 'Set2') +
  theme_classic() +
  theme(legend.position = 'none') +
  ggtitle('Class assignments before and after mutation\nGene B1, Positive regulation, log-sampled parameters\nWith basal paralog expression') +
  xlab('Mutated alleles') +
  ylab('Number of parameter sets')
ggsave(classes_sankey61_B1forfig, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/classes_sankey_B1formainfig.pdf'))


basic_class_assignment_all_forpie61 <- basic_class_assignment_all61 %>%
  group_by(mutated_alleles, product, class_assignment) %>%
  summarise(nSets = length(product))

classes_pies <- ggplot(basic_class_assignment_all_forpie61, aes(x="", y=nSets, fill=class_assignment)) +
  geom_bar(stat='identity', width=1, color='white') +
  coord_polar('y', start=0) +
  facet_grid(product~mutated_alleles) +
  theme_void() +
  ggtitle('classes of all distributions in v1.6.6.1\neach gene in each genotype')
ggsave(classes_pies, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/classes_pies.pdf'))

# get data into format: each row is a paramset, with cols: sampled LHS (numeric) and classes (factor) of each gene in each genotype

basic_class_assignment_all61_wide <- basic_class_assignment_all61 %>% 
  dplyr::select(version, paramset, mutated_alleles, product, class_assignment) %>%
  pivot_wider(names_from = c('mutated_alleles', 'product'), values_from = 'class_assignment')

classes_for_trees61 <- inner_join(basic_class_assignment_all61_wide, lhs_sets61, by = 'paramset')

temp_class_for_tree61 <- classes_for_trees61 %>%
  dplyr::select(colnames(lhs_sets61), `1_B1`) %>% ungroup() %>%
  dplyr::select(-paramset) %>%
  mutate(is_bimodal = ifelse(`1_B1` == 'bimodal', T, F)) %>%
  dplyr::select(-`1_B1`)

temp_class_for_tree61$is_bimodal <- as.factor(temp_class_for_tree61$is_bimodal)

temp.tree <- ctree(is_bimodal ~ basal_nitc_on_ratio + onbasalA1_off_ratio + A1_Aprime1_addon_ratio + 
                     A1_Aprime_prodon_ratio + r_prod_on + r_addon_byA1_B1 + r_onbasal_A1, 
                   data = temp_class_for_tree61,
                   control = ctree_control(alpha = 0.01))


# decision tree for B1 unimodal symm to unimodal symm vs unimodal symm to other


temp_class_for_tree61 <- classes_for_trees61 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  mutate(is_robust = ifelse(`0_B1` == `1_B1`, T, F)) %>%
  dplyr::select(colnames(lhs_sets61), is_robust) %>% ungroup() %>%
  dplyr::select(-paramset) 

temp_class_for_tree61$is_robust <- as.factor(temp_class_for_tree61$is_robust)

temp.tree <- ctree(is_robust ~ basal_nitc_on_ratio + onbasalA1_off_ratio + A1_Aprime1_addon_ratio + 
                     A1_Aprime_prodon_ratio + r_prod_on + r_addon_byA1_B1 + r_onbasal_A1, 
                   data = temp_class_for_tree61,
                   control = ctree_control(alpha = 0.01))

pdf(paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/isRobust_unimodalsymmetric_tree.pdf'), width = 20, height = 10)
plot(temp.tree)
dev.off()

# correlations

genos = 0:2
prods = c('A1', 'Anonsense1', 'Aprime1', 'B1')

pdf(paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/correlations_perGene_perGenotype.pdf'), width = 24, height = 32)
par(mfcol = c(4,3))
for (aci in 1:length(genos)) {
  for (proi in 1:length(prods)) {
    
    ac = genos[aci]
    pro = prods[proi]
    
    tempforcorr <- allstats_full61 %>% 
      filter(mutated_alleles == ac,
             product == pro) %>%
      ungroup() %>%
      inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>%
      dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'version', 'paramset', 'mutated_alleles', 'product'))
    
    corrplot(cor(tempforcorr), title = paste0(pro, ' in genotype ', as.character(ac)), mar = c(0,0,2,0))
    
  }
}
dev.off()

comparisons = unique(compared_stats61$compare)
for(comp in comparisons) {
  pdf(paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/correlations_',comp,'_perGene_perGenotype.pdf'), width = 12, height = 12)
  par(mfrow = c(2,2))
  for (proi in 1:length(prods)) {
    
    pro = prods[proi]
    
    tempforcorr <- compared_stats61 %>% 
      filter(compare == comp) %>% 
      dplyr::select(-c(`0`, `1`, `2`, 'mean_denom')) %>% 
      pivot_wider(names_from = stat, values_from = diff) %>% 
      filter(product == pro) %>%
      ungroup() %>%
      inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>%
      dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'version', 'paramset', 'product', 'compare'))
    
    corrplot(cor(tempforcorr), title = paste0(pro, ' compared ', comp), mar = c(0,0,2,0))
    
  }
  dev.off()
}


for (aci in 1:length(genos)) {
  
  for (proi in 1:length(prods)) {
    
    ac = genos[aci]
    pro = prods[proi]
    
    cat(paste0('Working on scatters for ', pro, ' in genotype ', ac, '\n'))
    
    tempforcorr <- allstats_full61 %>% 
      filter(mutated_alleles == ac,
             product == pro) %>%
      ungroup() %>%
      inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>%
      mutate(log10_mean_product1 = log10(mean_product+1)) %>%
      relocate(log10_mean_product1, .after = mean_product) %>%
      dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'mean_product')) %>%
      pivot_longer(names_to = 'statistic', values_to = 'stat_value', cols = log10_mean_product1:entropy) %>%
      pivot_longer(names_to = 'parameter', values_to = 'param_value', cols = basal_nitc_on_ratio:r_onbasal_Aprime_ratio)
    
    sc1 <- ggplot(tempforcorr, aes(log10(param_value), stat_value)) +
      geom_point(alpha = 0.1, stroke = 0) +
      facet_grid(statistic ~ parameter, scales = 'free') + 
      theme_classic() +
      xlab('log10(parameter value)') +
      ylab('statistic value') +
      ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles'))
    
    
    tempforcorr1 <- allstats_full61 %>% 
      filter(mutated_alleles == ac,
             product == pro,
             mean_product > 10) %>%
      ungroup() %>%
      inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>%
      mutate(log10_mean_product1 = log10(mean_product+1)) %>%
      relocate(log10_mean_product1, .after = mean_product) %>%
      dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'mean_product')) %>%
      pivot_longer(names_to = 'statistic', values_to = 'stat_value', cols = log10_mean_product1:entropy) %>%
      pivot_longer(names_to = 'parameter', values_to = 'param_value', cols = basal_nitc_on_ratio:r_onbasal_Aprime_ratio)
    ggsave(sc1, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatter_gene', pro, '_mutatedalleles', as.character(ac),'.pdf'), width = 24, height = 18)
    
    if(nrow(tempforcorr1)>0){
      sc1d <- ggplot(tempforcorr, aes(log10(param_value), stat_value)) +
        geom_point(alpha = 0.1, stroke = 0) +
        geom_density2d() +
        facet_grid(statistic ~ parameter, scales = 'free') + 
        theme_classic() +
        xlab('log10(parameter value)') +
        ylab('statistic value') +
        ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles'))
      ggsave(sc1d, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatterDensity_gene', pro, '_mutatedalleles', as.character(ac),'.pdf'), width = 24, height = 18)
      
      sc2 <- ggplot(tempforcorr1, aes(log10(param_value), stat_value)) +
        geom_point(alpha = 0.1, stroke = 0) +
        facet_grid(statistic ~ parameter, scales = 'free') + 
        theme_classic() +
        xlab('log10(parameter value)') +
        ylab('statistic value') +
        ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles\nMinimum mean expression = 10'))
      sc2d <- ggplot(tempforcorr1, aes(log10(param_value), stat_value)) +
        geom_point(alpha = 0.1, stroke = 0) +
        geom_density2d() +
        facet_grid(statistic ~ parameter, scales = 'free') + 
        theme_classic() +
        xlab('log10(parameter value)') +
        ylab('statistic value') +
        ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles\nMinimum mean expression = 10'))
      
      ggsave(sc2, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatter_minMean10_gene', pro, '_mutatedalleles', as.character(ac),'.pdf'), width = 24, height = 18)
      ggsave(sc2d, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatterDensity_minMean10_gene', pro, '_mutatedalleles', as.character(ac),'.pdf'), width = 24, height = 18)
    }
  }
}

# focus on bimodality coefficient vs paralog-relevant parameters
ac = 1
pro = 'B1'

cat(paste0('Working on scatters for ', pro, ' in genotype ', ac, '\n'))

tempforcorr <- allstats_full61 %>% 
  filter(mutated_alleles == ac,
         product == pro) %>%
  ungroup() %>%
  inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>%
  mutate(log10_mean_product1 = log10(mean_product+1)) %>%
  relocate(log10_mean_product1, .after = mean_product) %>%
  dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'mean_product')) %>%
  pivot_longer(names_to = 'statistic', values_to = 'stat_value', cols = log10_mean_product1:entropy) %>%
  pivot_longer(names_to = 'parameter', values_to = 'param_value', cols = basal_nitc_on_ratio:r_onbasal_Aprime_ratio) %>%
  filter(statistic == 'bimodality_coef',
         parameter %in% c('A1_Aprime1_addon_ratio', 'A1_Aprime_prodon_ratio', 'basal_nitc_on_ratio', 'r_onbasal_Aprime_ratio'))

sc1 <- ggplot(tempforcorr, aes(log10(param_value), stat_value)) +
  geom_point(alpha = 0.1, stroke = 0) +
  facet_grid(statistic ~ parameter, scales = 'free') + 
  theme_classic() +
  xlab('log10(parameter value)') +
  ylab('statistic value') +
  ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles'))


tempforcorr1 <- allstats_full61 %>% 
  filter(mutated_alleles == ac,
         product == pro,
         mean_product > 10) %>%
  ungroup() %>%
  inner_join(lhs_sets_full61, by = c('version', 'paramset')) %>%
  mutate(log10_mean_product1 = log10(mean_product+1)) %>%
  relocate(log10_mean_product1, .after = mean_product) %>%
  dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'mean_product')) %>%
  pivot_longer(names_to = 'statistic', values_to = 'stat_value', cols = log10_mean_product1:entropy) %>%
  pivot_longer(names_to = 'parameter', values_to = 'param_value', cols = basal_nitc_on_ratio:r_onbasal_Aprime_ratio) %>%
  filter(statistic == 'bimodality_coef',
         parameter %in% c('A1_Aprime1_addon_ratio', 'A1_Aprime_prodon_ratio', 'basal_nitc_on_ratio', 'r_onbasal_Aprime_ratio'))
ggsave(sc1, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatter_gene', pro, '_mutatedalleles', as.character(ac),'_B1hetfocus.pdf'), width = 16, height = 4)

if(nrow(tempforcorr1)>0){
  sc1d <- ggplot(tempforcorr, aes(log10(param_value), stat_value)) +
    geom_point(alpha = 0.1, stroke = 0) +
    geom_density2d() +
    facet_grid(statistic ~ parameter, scales = 'free') + 
    theme_classic() +
    xlab('log10(parameter value)') +
    ylab('statistic value') +
    ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles'))
  ggsave(sc1d, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatterDensity_gene', pro, '_mutatedalleles', as.character(ac),'_B1hetfocus.pdf'), width = 16, height = 4)
  
  sc2 <- ggplot(tempforcorr1, aes(log10(param_value), stat_value)) +
    geom_point(alpha = 0.1, stroke = 0) +
    facet_grid(statistic ~ parameter, scales = 'free') + 
    theme_classic() +
    xlab('log10(parameter value)') +
    ylab('statistic value') +
    ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles\nMinimum mean expression = 10'))
  sc2d <- ggplot(tempforcorr1, aes(log10(param_value), stat_value)) +
    geom_point(alpha = 0.1, stroke = 0) +
    geom_density2d() +
    facet_grid(statistic ~ parameter, scales = 'free') + 
    theme_classic() +
    xlab('log10(parameter value)') +
    ylab('statistic value') +
    ggtitle(paste0('Statistic vs parameter value\nGene: ', pro, ', Genotype: ', as.character(ac), ' mutated alleles\nMinimum mean expression = 10'))
  
  ggsave(sc2, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatter_minMean10_gene', pro, '_mutatedalleles', as.character(ac),'_B1hetfocus.pdf'), width = 16, height = 4)
  ggsave(sc2d, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatterDensity_minMean10_gene', pro, '_mutatedalleles', as.character(ac),'_B1hetfocus.pdf'), width = 16, height = 4)
}

cor(tempforcorr1 %>% 
      dplyr::filter(statistic == 'bimodality_coef', parameter == 'basal_nitc_on_ratio') %>% 
      dplyr::select(stat_value, param_value) %>% 
      mutate(log_param_val = log10(param_value)))

tempforcortest1 <- tempforcorr1 %>%
  pivot_wider(names_from = 'statistic', values_from = 'stat_value') %>%
  pivot_wider(names_from = 'parameter', values_from = 'param_value')

c_addon <- cor.test(tempforcortest1$bimodality_coef, log10(tempforcortest1$A1_Aprime1_addon_ratio))
c_absaddon <- cor.test(tempforcortest1$bimodality_coef, abs(log10(tempforcortest1$A1_Aprime1_addon_ratio)))
c_prodon <- cor.test(tempforcortest1$bimodality_coef, log10(tempforcortest1$A1_Aprime_prodon_ratio))
c_absprodon <- cor.test(tempforcortest1$bimodality_coef, abs(log10(tempforcortest1$A1_Aprime_prodon_ratio)))
c_nitc <- cor.test(tempforcortest1$bimodality_coef, log10(tempforcortest1$basal_nitc_on_ratio))
c_absnitc <- cor.test(tempforcortest1$bimodality_coef, abs(log10(tempforcortest1$basal_nitc_on_ratio)))
c_primon <- cor.test(tempforcortest1$bimodality_coef, log10(tempforcortest1$r_onbasal_Aprime_ratio))
c_absprimon <- cor.test(tempforcortest1$bimodality_coef, abs(log10(tempforcortest1$r_onbasal_Aprime_ratio)))

cortab <- tibble(
  statistic = c('A1_Aprime1_addon_ratio', 'abs(A1_Aprime1_addon_ratio)', 'A1_Aprime_prodon_ratio', 'abs(A1_Aprime_prodon_ratio)', 
                'basal_nitc_on_ratio', 'abs(basal_nitc_on_ratio)', 'r_onbasal_Aprime_ratio', 'abs(r_onbasal_Aprime_ratio)'),
  cor = c(c_addon$estimate, c_absaddon$estimate, c_prodon$estimate, c_absprodon$estimate, c_nitc$estimate, c_absnitc$estimate, c_primon$estimate, c_absprimon$estimate),
  p.val = c(c_addon$p.value, c_absaddon$p.value, c_prodon$p.value, c_absprodon$p.value, c_nitc$p.value, c_absnitc$p.value, c_primon$p.value, c_absprimon$p.value)
)
write.table(cortab, file = paste0(plotdir61, 'stats_class_assignment_check_v', as.character(anver),'/stats_scatterDensity_minMean10_gene', pro, '_mutatedalleles', as.character(ac),'_pearsonCor.txt'), quote = F, row.names = F)

