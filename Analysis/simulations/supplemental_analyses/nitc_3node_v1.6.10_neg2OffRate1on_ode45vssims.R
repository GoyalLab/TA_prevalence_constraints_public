# compares ode45-based steady-state analysis to simulations for 100 parameter sets

library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')

# edit as needed
datadir <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.10/samples/'
plotdir <- '~/code/grn_nitc/nitc_3node_v1.6.10/exploratory_analysis/'
setwd(datadir)

if(!dir.exists(plotdir)){
  dir.create(plotdir)
}

lhs_sets <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets %<>%
  mutate(paramset = 1:nrow(lhs_sets))

paramsets <- 1:100
allstats <- list()
allparams <- list()
for (paramset in paramsets){
  
  cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T)) #%>%
    # mutate(paramset = paramset,
    #        ssA_wt = 2*r_onbasal_A1*r_prodon_A1/(r_deg_A1*(r_onbasal_A1+r_off_A1)),
    #        HBA_wt = r_bind_byA1_B1*(ssA_wt^n_A1)/(k_A1^n_A1 + ssA_wt^n_A1),
    #        boundA_wt = HBA_wt/(HBA_wt+r_unbind_byA1_B1),
    #        onB_wt = r_bound_byA1_B1*boundA_wt/(r_bound_byA1_B1*boundA_wt + r_off_B1),
    #        ssB_wt = 2*r_prodon_B1*onB_wt/r_deg_B1)
  
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
  
  # dist_plot<-ggplot(species_sample, aes(abundance)) +
  #   geom_histogram() + 
  #   facet_grid(mutated_alleles~product) +
  #   ggtitle(paste0('Parameter set ', as.character(paramset))) +
  #   theme_classic()
  
  # ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_paramset_', as.character(paramset), '.pdf'))
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance))
  
  if(is.null(dim(allstats))) {
    allstats <- spstats
  } else {
    allstats %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams))) {
    allparams <- params
  } else {
    allparams %<>% bind_rows(params)
  }
  
}

#wtwt ode45
setwd(plotdir)
steadystate_ode45_wtwt <- as_tibble(read.csv('../steady_state_ODE45_wtwt_B3state_neg2OffRate1onODE.csv', header = T))
ssA_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 0) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtwt, by = 'paramset'), aes(A1, ss_A1)) +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nWT/WT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1') 

ssAn_plot <- ggplot(inner_join(allstats %>% 
                                 filter(mutated_alleles == 0) %>% 
                                 dplyr::select(paramset, product, mean_product) %>% 
                                 group_by(paramset) %>% 
                                 pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtwt, by = 'paramset'), aes(Anonsense1, ss_Anons1)) +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of Anons1 vs simulation\nWT/WT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Anonsense1') +
  ylab('ODE45 steady-state Anonsense1') 

ssAp_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 0) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtwt, by = 'paramset'), aes(Aprime1, ss_Aprim1)) +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of Aprim1 vs simulation\nWT/WT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Aprime1') +
  ylab('ODE45 steady-state Aprime1') 

ssB_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 0) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtwt, by = 'paramset'), aes(B1, ss_B1)) +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nWT/WT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_wtwt.pdf'), width = 10, height = 10)
ss_wt<-grid.arrange(ssA_plot, ssAn_plot, ssAp_plot, ssB_plot, ncol=2)
dev.off()

ssB_plot_n <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 2) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset') %>% inner_join(allparams %>% mutate(paramset = 1:100), by = 'paramset'), aes(B1, ss_B1)) +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_point() +
  # geom_text_repel(aes(label = as.character(round(n_B1,2)))) +
  theme_bw() +
  xlim(c(0,500)) + ylim(c(0,500)) +
  ggtitle('Steady-state approximation of B1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

## het ode45 from matlab
steadystate_ode45_wtmut <- as_tibble(read.csv('../steady_state_ODE45_wtmut_B3state_neg2OffRate1onODE.csv', header = T))
ssA_wm_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 1) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(A1, ss_A1)) +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1')
ssAnons_wm_plot <- ggplot(inner_join(allstats %>% 
                                       filter(mutated_alleles == 1) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(Anonsense1, ss_Anons1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of Anonsense1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Anonsense1') +
  ylab('ODE45 steady-state Anonsense1')
ssAprim_wm_plot <- ggplot(inner_join(allstats %>% 
                                       filter(mutated_alleles == 1) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(Aprime1, ss_Aprim1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of Aprime1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Aprime1') +
  ylab('ODE45 steady-state Aprime1')
ssB_wm_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 1) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_wtmut.pdf'), width = 10, height = 10)
ss_wt<-grid.arrange(ssA_wm_plot, ssAnons_wm_plot, ssAprim_wm_plot, ssB_wm_plot, ncol=2)
dev.off()



steadystate_ode45_mutmut <- as_tibble(read.csv('../steady_state_ODE45_mutmut_B3state_neg2OffRate1onODE.csv', header = T))
ssA_mm_plot <- ggplot(inner_join(allstats %>% 
                                   filter(mutated_alleles == 2) %>% 
                                   dplyr::select(paramset, product, mean_product) %>% 
                                   group_by(paramset) %>% 
                                   pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(A1, ss_A1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1')
ssAnons_mm_plot <- ggplot(inner_join(allstats %>% 
                                       filter(mutated_alleles == 2) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(Anonsense1, ss_Anons1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of Anonsense1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Anonsense1') +
  ylab('ODE45 steady-state Anonsense1')
ssAprim_mm_plot <- ggplot(inner_join(allstats %>% 
                                       filter(mutated_alleles == 2) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(Aprime1, ss_Aprim1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of Aprime1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Aprime1') +
  ylab('ODE45 steady-state Aprime1')
ssB_mm_plot <- ggplot(inner_join(allstats %>% 
                                   filter(mutated_alleles == 2) %>% 
                                   dplyr::select(paramset, product, mean_product) %>% 
                                   group_by(paramset) %>% 
                                   pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_mutmut.pdf'), width = 10, height = 10)
ss_wt<-grid.arrange(ssA_mm_plot, ssAprim_mm_plot, ssAnons_mm_plot, ssB_mm_plot, ncol=2)
dev.off()

ssB_mm_plot_params <- ggplot(inner_join(allstats %>% 
                                          filter(mutated_alleles == 2) %>% 
                                          dplyr::select(paramset, product, mean_product) %>% 
                                          group_by(paramset) %>% 
                                          pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset') %>%
                               inner_join(allparams %>% mutate(paramset = 1:100), by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'blue', linetype = 2) +
  geom_text_repel(aes(label = ifelse(abs(log(ss_B1/B1))>0.5 & B1>20, as.character(round(n_B1, 2)), ''))) +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

off_diag_ssB_mutmut <- inner_join(allstats %>% 
                                    filter(mutated_alleles == 2) %>% 
                                    dplyr::select(paramset, product, mean_product) %>% 
                                    group_by(paramset) %>% 
                                    pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset') %>%
  filter(abs(log(ss_B1/B1))>0.5, B1>20)

lhs_sets_offdiag <- lhs_sets %>%
  mutate(off_diag_ssB_mutmut = paramset %in% off_diag_ssB_mutmut$paramset)

logit1<-glm(off_diag_ssB_mutmut ~ ., data = lhs_sets_offdiag, family = 'binomial')


lhs_sets_offdiag2 <- lhs_sets %>%
  mutate(off_diag_ssB_mutmut = (paramset %in% off_diag_ssB_mutmut$paramset & paramset != 99))

logit2<-glm(off_diag_ssB_mutmut ~ ., data = lhs_sets_offdiag2, family = 'binomial')

# make traces for off_diag paramsets
off_diag_ssB_mutmut
on_diag_ssB_mutmut <- lhs_sets %>%
  inner_join(allstats %>% filter(product == 'B1', mutated_alleles == 2), by = 'paramset') %>%
  filter(!(paramset %in% off_diag_ssB_mutmut$paramset), mean_product > 24, mean_product < 107)

offdiagdir <- paste0(plotdir, 'off_diagonal_B1_mutmut/')
ondiagdir <- paste0(plotdir, 'on_diagonal_B1_mutmut/')

if(!dir.exists(offdiagdir)){ dir.create(offdiagdir)}
if(!dir.exists(ondiagdir)){ dir.create(ondiagdir)}
orig_color = 'black'
nons_color = 'firebrick2'
para_color = 'gray'
targ_color = 'dodgerblue2'

maxtime = 210400

setwd(datadir)
for (pset in off_diag_ssB_mutmut$paramset) {
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(pset),'.csv'), header = T))
  nr = nrow(species)
  species %<>%
    mutate(time = 1:nrow(species),
           paramset = paramset) %>%
    filter(time > 200400, time < maxtime)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.3) +
    geom_line(data = species, aes(time, B1), color = targ_color) +
    geom_hline(data = off_diag_ssB_mutmut %>% ungroup() %>% filter(paramset == pset), aes(yintercept=B1), color = 'blue', linetype = 1) +
    geom_hline(data = off_diag_ssB_mutmut %>% ungroup() %>% filter(paramset == pset), aes(yintercept=ss_B1), color = 'blue', linetype = 2) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, parameter set ', as.character(pset),'\nSolid line = sampled mean, Dashed line = steady state'))
  ggsave(spec_plot, file=paste0(offdiagdir, 'species_trace_mutmut_', as.character(pset), '.pdf'), width = 10, height = 3)
  
}

on_diag_ssB_mutmut_stats<-inner_join(allstats %>% ungroup() %>%
                                       filter(mutated_alleles == 2) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset') %>%
  filter(paramset %in% on_diag_ssB_mutmut$paramset)
for (pset in on_diag_ssB_mutmut$paramset) {
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(pset),'.csv'), header = T))
  nr = nrow(species)
  species %<>%
    mutate(time = 1:nrow(species),
           paramset = paramset) %>%
    filter(time > 200400, time < maxtime)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.3) +
    geom_line(data = species, aes(time, B1), color = targ_color) +
    geom_hline(data = on_diag_ssB_mutmut_stats %>% ungroup() %>% filter(paramset == pset), aes(yintercept=B1), color = 'blue', linetype = 1) +
    geom_hline(data = on_diag_ssB_mutmut_stats %>% ungroup() %>% filter(paramset == pset), aes(yintercept=ss_B1), color = 'blue', linetype = 2) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, parameter set ', as.character(pset),'\nSolid line = sampled mean, Dashed line = steady state'))
  ggsave(spec_plot, file=paste0(ondiagdir, 'species_trace_mutmut_', as.character(pset), '.pdf'), width = 10, height = 3)
  
}  


# ode45 with consecutive starting conditions (like in sims)

#wtwt ode45
steadystate_ode45_wtwt <- as_tibble(read.csv('../steady_state_ODE45_wtwt3.csv', header = T))
ssA_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 0) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtwt, by = 'paramset'), aes(A1, ss_A1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nWT/WT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1')

ssB_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 0) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtwt, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nWT/WT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_wtwt_AB_consecutiveIC.pdf'), width = 10, height = 5)
ss_wt<-grid.arrange(ssA_plot, ssB_plot, ncol=2)
dev.off()

## het ode45 from matlab
steadystate_ode45_wtmut <- as_tibble(read.csv('../steady_state_ODE45_wtmut3.csv', header = T))
ssA_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 1) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(A1, ss_A1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1')

ssB_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 1) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_wtmut_AB_consecutiveIC.pdf'), width = 10, height = 5)
ss_wt<-grid.arrange(ssA_plot, ssB_plot, ncol=2)
dev.off()



steadystate_ode45_mutmut <- as_tibble(read.csv('../steady_state_ODE45_mutmut3.csv', header = T))
ssA_mm_plot <- ggplot(inner_join(allstats %>% 
                                   filter(mutated_alleles == 2) %>% 
                                   dplyr::select(paramset, product, mean_product) %>% 
                                   group_by(paramset) %>% 
                                   pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(A1, ss_A1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1')
ssAnons_mm_plot <- ggplot(inner_join(allstats %>% 
                                       filter(mutated_alleles == 2) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(Anonsense1, ss_Anons1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of Anonsense1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Anonsense1') +
  ylab('ODE45 steady-state Anonsense1')
ssAprim_mm_plot <- ggplot(inner_join(allstats %>% 
                                       filter(mutated_alleles == 2) %>% 
                                       dplyr::select(paramset, product, mean_product) %>% 
                                       group_by(paramset) %>% 
                                       pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(Aprime1, ss_Aprim1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of Aprime1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean Aprime1') +
  ylab('ODE45 steady-state Aprime1')
ssB_mm_plot <- ggplot(inner_join(allstats %>% 
                                   filter(mutated_alleles == 2) %>% 
                                   dplyr::select(paramset, product, mean_product) %>% 
                                   group_by(paramset) %>% 
                                   pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_mutmut_AB_consecutiveIC.pdf'), width = 10, height = 10)
ss_wt<-grid.arrange(ssA_mm_plot, ssAprim_mm_plot, ssAnons_mm_plot, ssB_mm_plot, ncol=2)
dev.off()

ssB_mm_plot_params <- ggplot(inner_join(allstats %>% 
                                          filter(mutated_alleles == 2) %>% 
                                          dplyr::select(paramset, product, mean_product) %>% 
                                          group_by(paramset) %>% 
                                          pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_mutmut, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  geom_text_repel(aes(label = ifelse(abs(log(ss_B1/B1))>0.5 & B1>20, paramset, ''))) +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nMUT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')
