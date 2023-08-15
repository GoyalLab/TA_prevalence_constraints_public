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
datadir53 <- '/Volumes/IAMYG2/grn_nitc_data/v1.6.5.3/samples/'
plotdir53 <- '/Volumes/IAMYG2/grn_nitc_data/v1.6.5.3/exploratory_analysis/'
setwd(datadir53)

if(!dir.exists(plotdir53)){
  dir.create(plotdir53)
}
paramsets53 <- 1:100
lhs_sets53 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets53 %<>%
  mutate(paramset = 1:nrow(lhs_sets53))

# calculate stats
allstats53 <- list()
allparams53 <- list()
for (paramset in paramsets53){
  
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
    ggsave(dist_plot, file = paste0(plotdir53, 'distributions_q300_v1.6.5.3_paramset_', as.character(paramset), '.pdf'))
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
  
  if(is.null(dim(allstats53))) {
    allstats53 <- spstats
  } else {
    allstats53 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams53))) {
    allparams53 <- params
  } else {
    allparams53 %<>% bind_rows(params)
  }
  
}


# temp: collate all data
setwd(datadir53)
all_species_q300_53 <- list()
for (paramset in paramsets53){
  
  if(paramset %% 100 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
  }
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset,
           version = '1.6.5.3') %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, version, time, mutated_alleles)
  
  
  if(is.null(dim(all_species_q300_53))) {
    all_species_q300_53 <- species_sample
  } else {
    all_species_q300_53 %<>% bind_rows(species_sample)
  }
  
}

write.csv(all_species_q300_53, file = paste0('/Volumes/IAMYG2/grn_nitc_data/v1.6.5.3/all_species_q300.csv'))


# summary stats draft1
allstats_full53 <- allstats53 %>% mutate(version = '1.6.5.3')
allparams_full53 <- allparams53 %>% mutate(paramset = 1:100, version = '1.6.5.3')
lhs_sets_full53 <- lhs_sets53 %>% mutate(version = '1.6.5.3')

pseud = 0.01

allstats_full53 %<>% mutate(skewness = ifelse(is.na(skewness), 0, skewness))

compared_stats53 <- allstats_full53 %>% 
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
  inner_join(allstats_full53 %>% 
               dplyr::select(mutated_alleles, product, paramset, version, mean_product) %>%
               pivot_wider(names_from = mutated_alleles, values_from = mean_product), by = c('product', 'paramset', 'version')) %>%
  mutate(mean_denom = case_when(
    compare %in% c('lfc10', 'delta10', 'lfc20', 'delta20') ~ `0`,
    compare %in% c('lfc21', 'delta21') ~ `1`))

# filter to Hill n < 5 - now already done for allstats_full1
#compared_stats %<>% ungroup() %>% inner_join(lhs_sets_fall %>% filter(Hill_coefficient_n < 5) %>% dplyr::select(version, paramset), by = c('version','paramset'))

# plot wt/mut summary stats against mean expression

unistats53<-unique(compared_stats53$stat)

for (st in unistats53) {
  
  pvs1 <- ggplot(allstats_full53 %>% inner_join(lhs_sets_full53, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n<5)) + 
    geom_point(aes(mean_product, eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = 10), linetype = 2) +
    # geom_text(aes(mean_product, eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic() +
    ggtitle(paste0(st, ' vs mean'))
  
  pvs2 <- ggplot(allstats_full53 %>% inner_join(lhs_sets_full53, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n<5)) + 
    geom_point(aes(log(mean_product), eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = log(10)), linetype = 2) +
    # geom_text(aes(log(mean_product), eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic()  +
    ggtitle(paste0(st, ' vs mean'))
  
  pvs3 <- ggplot(allstats_full53 %>% 
                   inner_join(lhs_sets_full53, by = c('version', 'paramset')) %>% 
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
  
  ggsave(pvs1, file = paste0(plotdir53, 'PerGenotype_', st, '_vs_mean.pdf'), width = 16, height = 8) 
  ggsave(pvs2, file = paste0(plotdir53, 'PerGenotype_', st, '_vs_logmean.pdf'), width = 16, height = 8) 
  ggsave(pvs3, file = paste0(plotdir53, 'WTMUT_', st, '_vs_logmean.pdf'), width = 4, height = 4) 
}


# LOESS
loess_fitted_allstats_all53 <- allstats_full53
for (stat in unistats53[unistats53 != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  statdat <- list()
  statdat1 <- list()
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    # for (ma in 0:2) {
    
    tempdat <- allstats_full53 %>% 
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
    
    # ggsave(lplot1, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.2and5.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_v1.6.2only.pdf'), width = 5, height = 5)
    # ggsave(lplot2, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.2and5.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_log_v1.6.2only.pdf'), width = 5, height = 5)
    
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
  
  
  loess_fitted_allstats_all53 %<>% left_join(as_tibble(statdat) %>% dplyr::select(-mean_product), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  # 
  # cat('sliding window normalizing...\n') # do this per-gene over all genotypes, rather than all genes over all genotypes...
  # statdat1 <- sliding_window_normalize(as_tibble(statdat) %>% filter(mean_product>10), 'mean_product', paste0(stat,'_residual'), 50)
  # 
  
  loess_fitted_allstats_all53 %<>% left_join(statdat1 %>% dplyr::select(-c('mean_product', paste0(stat,'_residual'), paste0(stat,'_fitted'))), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  
  cat(paste0('Done with ', stat, '\n'))
  
}

write.csv(loess_fitted_allstats_all53, file = paste0(plotdir53, 'loess_fitted_allstats_all_snwRadius100.csv'), quote = F, row.names = F)
loess_fitted_allstats_all53 <- as_tibble(read.csv(paste0(plotdir53, 'loess_fitted_allstats_all_snwRadius100.csv')))
# load as needed
# loess_fitted_allstats_all53 <- as_tibble(read.csv(paste0(plotdir53, 'loess_fitted_allstats_all_snwRadius100.csv'), header=T, stringsAsFactors = F))

for (stat in unistats53[unistats53 != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  # statdat <- list()
  
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    
    statdatA = loess_fitted_allstats_all53 %>%
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
    
    pdf(paste0(paste0(plotdir53, 'LOESSplots_', stat, 'vsLogMean_', gene,'_v1.6.2and5.pdf')), width = 10, height = 7)
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

basic_class_assignment_all53 <- loess_fitted_allstats_all53 %>%
  mutate(class_assignment = case_when(
    mean_product < 10 ~ 'low-average',
    bimodality_coef_residual > bimfilt & bimodality_coef > 0.555 ~ 'bimodal',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & abs(skewness) < 1 ~ 'unimodal symmetric',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness >= 1 ~ 'right-skewed unimodal',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness <= -1 ~ 'left-skewed unimodal'
    
  )) 

basic_class_assignment_all_forSankey53 <- basic_class_assignment_all53 %>%
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
basic_class_assignment_all_forSankey_forsamples53 <- basic_class_assignment_all53 %>%
  dplyr::select(version, paramset, product, mutated_alleles, class_assignment)  

classes53 = unique(basic_class_assignment_all53$class_assignment)

if(!dir.exists(paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver)))) {
  dir.create(paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver)))
}

classes_sankey <- ggplot(basic_class_assignment_all_forSankey53, aes(x = mutated_alleles, y=Freq,
                                                                      stratum = class_assignment, alluvium = alluvID, fill = class_assignment, label = class_assignment)) +
  facet_grid(product~.) +
  geom_flow() +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = 'stratum', size = 3) + 
  theme(legend.position = 'none')
ggsave(classes_sankey, file = paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/classes_sankey.pdf'))

basic_class_assignment_all_forpie53 <- basic_class_assignment_all53 %>%
  group_by(mutated_alleles, product, class_assignment) %>%
  summarise(nSets = length(product))

classes_pies <- ggplot(basic_class_assignment_all_forpie53, aes(x="", y=nSets, fill=class_assignment)) +
  geom_bar(stat='identity', width=1, color='white') +
  coord_polar('y', start=0) +
  facet_grid(product~mutated_alleles) +
  theme_void() +
  ggtitle('classes of all distributions in v1.6.5.3\neach gene in each genotype')
ggsave(classes_pies, file = paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/classes_pies.pdf'))

# get data into format: each row is a paramset, with cols: sampled LHS (numeric) and classes (factor) of each gene in each genotype

basic_class_assignment_all53_wide <- basic_class_assignment_all53 %>% 
  dplyr::select(version, paramset, mutated_alleles, product, class_assignment) %>%
  pivot_wider(names_from = c('mutated_alleles', 'product'), values_from = 'class_assignment')

classes_for_trees53 <- inner_join(basic_class_assignment_all53_wide, lhs_sets53, by = 'paramset')

temp_class_for_tree53 <- classes_for_trees53 %>%
  dplyr::select(colnames(lhs_sets53), `1_B1`) %>% ungroup() %>%
  dplyr::select(-paramset) %>%
  mutate(is_bimodal = ifelse(`1_B1` == 'bimodal', T, F)) %>%
  dplyr::select(-`1_B1`)

temp_class_for_tree53$is_bimodal <- as.factor(temp_class_for_tree53$is_bimodal)

temp.tree <- ctree(is_bimodal ~ basal_nitc_on_ratio + onbasalA1_off_ratio + A1_Aprime1_addon_ratio + 
                     A1_Aprime_prodon_ratio + r_prod_on + r_addon_byA1_B1 + r_onbasal_A1, 
                   data = temp_class_for_tree53,
                   control = ctree_control(alpha = 0.01))


# decision tree for B1 unimodal symm to unimodal symm vs unimodal symm to other


temp_class_for_tree53 <- classes_for_trees53 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  mutate(is_robust = ifelse(`0_B1` == `1_B1`, T, F)) %>%
  dplyr::select(colnames(lhs_sets53), is_robust) %>% ungroup() %>%
  dplyr::select(-paramset) 

temp_class_for_tree53$is_robust <- as.factor(temp_class_for_tree53$is_robust)

temp.tree <- ctree(is_robust ~ basal_nitc_on_ratio + onbasalA1_off_ratio + A1_Aprime1_addon_ratio + 
                     A1_Aprime_prodon_ratio + r_prod_on + r_addon_byA1_B1 + r_onbasal_A1, 
                   data = temp_class_for_tree53,
                   control = ctree_control(alpha = 0.01))

pdf(paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/isRobust_unimodalsymmetric_tree.pdf'), width = 20, height = 10)
plot(temp.tree)
dev.off()

classes_for_trees53 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  inner_join(compared_stats53 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = 'paramset') %>%
  group_by(`1_B1`) %>%
  summarise(length(`1_B1`))

resampled_US_LFC53 <- classes_for_trees53 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  inner_join(compared_stats53 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = 'paramset') %>%
  filter(`1_B1` == 'unimodal symmetric') %>%
  ungroup() %>%
  summarise(meanLFC = mean(diff), 
            semLFC = sd(diff)/sqrt(length(diff))) %>%
  mutate(nodeID = 'resampled node 15+16') %>%
  bind_rows(unimodal_symmetric_robust_tree52 %>%
              inner_join(compared_stats52 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = 'paramset') %>%
              filter(`1_B1` == 'unimodal symmetric') %>%
              group_by(nodeID) %>%
              summarise(meanLFC = mean(diff),
                        semLFC = sd(diff)/sqrt(length(diff))))

lfc_aftermutation_barplot_52and53 <- ggplot() +
  geom_bar(data = resampled_US_LFC53, aes(nodeID, meanLFC), stat = 'identity') +
  geom_errorbar(data = resampled_US_LFC53, aes(x = nodeID, ymin = meanLFC - semLFC, ymax = meanLFC + semLFC), width = 0.25) +
  theme_classic() +
  ylab('mean LFC after mutation') +
  xlab('branch grouping') +
  ggtitle('Mean log fold change after mutation\nUnimodal symmetric robust parameter sets only') +
  theme(axis.text.x = element_text(angle = 45))

ggsave(lfc_aftermutation_barplot_52and53, file = paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/lfc_aftermutation_barplot_52and53.pdf'))

set.seed(48332)
classes_for_trees52_samp<-classes_for_trees52[sample(1:10000, 100),] 

n53us <- sum(classes_for_trees53$`0_B1` == 'unimodal symmetric')
n52us <- sum(classes_for_trees52_samp$`0_B1` == 'unimodal symmetric')
n5253us <- tibble(
  version = c('1.6.5.2', '1.6.5.3'),
  n_total = c(n52us, n53us)
)

unimodal_robust_compare5352 <- classes_for_trees53 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  dplyr::select(version, paramset, `1_B1`) %>%
  inner_join(compared_stats53 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = c('version', 'paramset')) %>%
  mutate(isRobust = (`1_B1` == 'unimodal symmetric') & (abs(diff)<0.25)) %>%
  bind_rows(classes_for_trees52_samp %>%
              filter(`0_B1` == 'unimodal symmetric') %>%
              dplyr::select(version, paramset, `1_B1`) %>%
              inner_join(compared_stats52 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = c('version', 'paramset')) %>%
              mutate(isRobust = (`1_B1` == 'unimodal symmetric') & (abs(diff)<0.25))) %>%
  group_by(version, isRobust) %>%
  summarise(n_each = length(isRobust)) %>%
    inner_join(n5253us) %>%
  mutate(fraction_each = n_each/n_total)

unimodal_robust_compare5352_bars <- ggplot(unimodal_robust_compare5352, aes(version, fraction_each, fill=isRobust)) +
  geom_bar(stat = 'identity') +
  theme_classic() 

ggsave(unimodal_robust_compare5352_bars, file = paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/unimodal_robust_withmean_compare5352_bars.svg'))

unimodal_robust_compare5352_015 <- classes_for_trees53 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  dplyr::select(version, paramset, `1_B1`) %>%
  inner_join(compared_stats53 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = c('version', 'paramset')) %>%
  mutate(isRobust = (`1_B1` == 'unimodal symmetric') & (abs(diff)<0.15)) %>%
  bind_rows(classes_for_trees52_samp %>%
              filter(`0_B1` == 'unimodal symmetric') %>%
              dplyr::select(version, paramset, `1_B1`) %>%
              inner_join(compared_stats52 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = c('version', 'paramset')) %>%
              mutate(isRobust = (`1_B1` == 'unimodal symmetric') & (abs(diff)<0.15))) %>%
  group_by(version, isRobust) %>%
  summarise(n_each = length(isRobust)) %>%
  inner_join(n5253us) %>%
  mutate(fraction_each = n_each/n_total)

unimodal_robust_compare5352_bars_015 <- ggplot(unimodal_robust_compare5352_015, aes(version, fraction_each, fill=isRobust)) +
  geom_bar(stat = 'identity') +
  theme_classic() 

ggsave(unimodal_robust_compare5352_bars_015, file = paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/unimodal_robust_withmean_compare5352_bars_015.svg'))

unimodal_robust_compare5352_035 <- classes_for_trees53 %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  dplyr::select(version, paramset, `1_B1`) %>%
  inner_join(compared_stats53 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = c('version', 'paramset')) %>%
  mutate(isRobust = (`1_B1` == 'unimodal symmetric') & (abs(diff)<0.35)) %>%
  bind_rows(classes_for_trees52_samp %>%
              filter(`0_B1` == 'unimodal symmetric') %>%
              dplyr::select(version, paramset, `1_B1`) %>%
              inner_join(compared_stats52 %>% filter(product == 'B1', stat == 'mean_product', compare == 'lfc10'), by = c('version', 'paramset')) %>%
              mutate(isRobust = (`1_B1` == 'unimodal symmetric') & (abs(diff)<0.35))) %>%
  group_by(version, isRobust) %>%
  summarise(n_each = length(isRobust)) %>%
  inner_join(n5253us) %>%
  mutate(fraction_each = n_each/n_total)

unimodal_robust_compare5352_bars_035 <- ggplot(unimodal_robust_compare5352_035, aes(version, fraction_each, fill=isRobust)) +
  geom_bar(stat = 'identity') +
  theme_classic() 

ggsave(unimodal_robust_compare5352_bars_035, file = paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/unimodal_robust_withmean_compare5352_bars_035.svg'))


#correlations

genos = 0:2
prods = c('A1', 'Anonsense1', 'Aprime1', 'B1')

pdf(paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/correlations_perGene_perGenotype.pdf'), width = 24, height = 32)
par(mfcol = c(4,3))
for (aci in 1:length(genos)) {
  for (proi in 1:length(prods)) {
    
    ac = genos[aci]
    pro = prods[proi]
    
    tempforcorr <- allstats_full53 %>% 
      filter(mutated_alleles == ac,
             product == pro) %>%
      ungroup() %>%
      inner_join(lhs_sets_full53, by = c('version', 'paramset')) %>%
      dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'version', 'paramset', 'mutated_alleles', 'product'))
    
    corrplot(cor(tempforcorr), title = paste0(pro, ' in genotype ', as.character(ac)), mar = c(0,0,2,0))
    
  }
}
dev.off()

comparisons = unique(compared_stats53$compare)
for(comp in comparisons) {
  pdf(paste0(plotdir53, 'stats_class_assignment_check_v', as.character(anver),'/correlations_',comp,'_perGene_perGenotype.pdf'), width = 12, height = 12)
  par(mfrow = c(2,2))
  for (proi in 1:length(prods)) {
    
    pro = prods[proi]
    
    tempforcorr <- compared_stats53 %>% 
      filter(compare == comp) %>% 
      dplyr::select(-c(`0`, `1`, `2`, 'mean_denom')) %>% 
      pivot_wider(names_from = stat, values_from = diff) %>% 
      filter(product == pro) %>%
      ungroup() %>%
      inner_join(lhs_sets_full53, by = c('version', 'paramset')) %>%
      dplyr::select(-c('sd_product', 'fano_product', 'entropy95', 'entropy90', 'version', 'paramset', 'product', 'compare'))
    
    corrplot(cor(tempforcorr), title = paste0(pro, ' compared ', comp), mar = c(0,0,2,0))
    
  }
  dev.off()
}

