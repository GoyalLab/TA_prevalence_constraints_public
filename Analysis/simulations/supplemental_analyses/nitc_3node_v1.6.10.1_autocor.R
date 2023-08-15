
library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')

# creates autocorrelation plots for a random sample of parameter sets at q100 timestep sample frequency


# edit as needed
datadir <- '/Volumes/IAMYG2/grn_nitc_data/v1.6.10.1/fullTraces/'
plotdir <- '/Volumes/IAMYG2/grn_nitc_data/v1.6.10.1/exploratory_analysis/autocor/'
setwd(datadir)

if(!dir.exists(plotdir)){
  dir.create(plotdir)
}

paramsets <- 1:3
for (paramset in paramsets){
  
  cat(paste0('Working on ', as.character(paramset), '\n'))
  
  # params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T)) #%>%
  # mutate(paramset = paramset,
  #        ssA_wt = 2*r_onbasal_A1*r_prodon_A1/(r_deg_A1*(r_onbasal_A1+r_off_A1)),
  #        HBA_wt = r_bind_byA1_B1*(ssA_wt^n_A1)/(k_A1^n_A1 + ssA_wt^n_A1),
  #        boundA_wt = HBA_wt/(HBA_wt+r_unbind_byA1_B1),
  #        onB_wt = r_bound_byA1_B1*boundA_wt/(r_bound_byA1_B1*boundA_wt + r_off_B1),
  #        ssB_wt = 2*r_prodon_B1*onB_wt/r_deg_B1)
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  
  steady_state_sample_wtwt <- species %>%
    filter(time > 400, time < 100000, (time + 1) %% 100 == 0) 
  
  steady_state_sample_wtmut <- species %>%
    filter(time > 100400, time < 200000, (time + 1) %% 100 == 0)
  
  steady_state_sample_mutmut <- species %>%
    filter(time > 200400, (time + 1) %% 100 == 0)
  
  autocor_wws<-acf(steady_state_sample_wtwt$B1)
  autocor_wms<-acf(steady_state_sample_wtmut$B1)
  autocor_mms<-acf(steady_state_sample_mutmut$B1)
  
  
  pdf(paste0(plotdir, 'autocor_steadstatesamples_q100_paramset', as.character(paramset), '.pdf'), height = 5, width = 15)
  par(mfrow=c(1,3))
  plot(autocor_wws, main = paste0('Autocorrelation per 100 time steps, WT/WT\nParameter set ', as.character(paramset)))
  plot(autocor_wms, main = paste0('Autocorrelation per 100 time steps, WT/MUT\nParameter set ', as.character(paramset)))
  plot(autocor_mms, main = paste0('Autocorrelation per 100 time steps, MUT/MUT\nParameter set ', as.character(paramset)))
  dev.off()
  
}


