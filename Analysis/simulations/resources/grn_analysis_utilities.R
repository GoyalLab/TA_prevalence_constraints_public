## Utility functions for analysis of GRN output

library(e1071)

# default color palette for classifier
# dist_pal_default matches colors from the colorblind-friendly Rcolorbrewer palette Set2 to classes in analysis version 5 classes

# calculate_bimodality calculates Sarle's bimodality coefficient for a finite sample (https://en.wikipedia.org/wiki/Multimodal_distribution#Bimodality_coefficient).
# Note that a unimodal distribution will tend toward 0, a uniform distribution will give a value of ~0.555, and a Bernoulli distribution will give a value of 1.
# Accepts:
# - x: a numeric vector.
#
# Returns:
# - bc: a numeric value of Sarle's bimodality coefficient, between 0 and 1. Will return 0 for spike-and-slab dist (e.g., when all values 0)
calculate_bimodality <- function(x) {
  
  nx <- length(x)
  
  if(sum(x == mean(x)) == nx) {
    bc = 0
  } else {
    
    sx <- skewness(x)
    kx <- kurtosis(x) # excess kurtosis = m_4 / m_2^2 - 3
    
    bc <- (sx^2 + 1) / (kx + (3*(nx-1)^2)/((nx-2)*(nx-3)))
  }
  return(bc)
  
}

# sliding_window_normalize performs sliding-window normalization to a 0-1 scale for a dependent variable along an independent variable axis (e.g., 0-1 scale normalize CV along mean axis)
# Accepts:
# - dat: a tibble with columns including the indep and depen column names
# - indep: a character, the colname of the independent variable (e.g., mean)
# - depen: a character, the colname of the dependent variable (e.g., CV)
# - radius: an integer, the sliding window radius (number of values upstream and number of values downstream to consider. E.g., radius = 25 gives a window of size 51 including the reference value)
#
# Returns:
# - dat_mod: a tibble, dat with new columns 'depen'_swn
sliding_window_normalize <- function(dat, indep, depen, radius) {
  
  nr <- nrow(dat)
  
  dat %<>% ungroup() %>% mutate(tID = 1:nr)
  
  tempdat <- dat %>% ungroup() %>% dplyr::select(as.symbol(eval(indep)), as.symbol(eval(depen)), tID)
  
  tempdat <- tempdat[order(tempdat$mean_product),] %>% as.data.frame()
  
  dat_mod <- list()
  for (i in 1:nr) {
    
    tind <- tempdat[i,indep]
    tdep <- tempdat[i,depen]
    tid <- tempdat[i,'tID']
    
    td1 <- tempdat[max(1, (i-radius)):min(nr, (i+radius)),]
    
    maxdepen = max(td1[,depen]) + 0.00001 # prevent divide by zero
    mindepen = min(td1[,depen])
    
    tdep_swn = (tdep-mindepen)/(maxdepen-mindepen)
    
    tdat_swn <- tibble(
      X1 = tind[1],
      X2 = tdep,
      X3 = tdep_swn,
      X4 = tid
    )
    
    if(is.null(dim(dat_mod))){
      dat_mod <- tdat_swn
    } else {
      dat_mod %<>% bind_rows(tdat_swn)
    }
    
  }
  
  colnames(dat_mod) <- c(indep, depen, paste0(depen, '_swn'), 'tID')
  
  dat_mod %<>% inner_join(dat, by = c('tID', eval(depen), eval(indep))) %>% dplyr::select(-tID)
  
  return(dat_mod)

}
  

# plot_traces plots example traces for data from a simulation given timepoint bounds.
# Accepts:
# - dat: a tibble with columns: time, A1, Anonsense1, Aprime1, B1, paramset
# - start: numeric, starting time
# - end: numeric, ending time
# - genostring: string, describing genotype, e.g., 'WT/WT'
# 
# Returns:
# - p1: ggplot with traces only
plot_traces <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.5) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.5) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

plot_traces_ver <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.5) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.5) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, version ', as.character(ver), ', parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

plot_traces_ver_Bfocus <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.2) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.9) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, version ', as.character(ver), ', parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

plot_traces_ver_nonons <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.5) +
    # geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.5) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.5) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, version ', as.character(ver), ', parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

plot_traces_ver_nonons_Bfocus <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.2) +
    # geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.2) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.9) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, version ', as.character(ver), ', parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

##
# tree - a ctree object
# paramspace - a tibble with columns 'parameter', 'minimum', 'maximum'
##
# return leafspace: a tibble defining parameter subspaces with columns 'leafID', 'parameter', 'lower', 'upper', 'range' (range calculated at end)
parse_partykit <- function(tree, paramspace) {
  
  noderules <- partykit:::.list.rules.party(temp.tree)
  
  leafspace <- tibble(
    leafID = numeric(),
    parameter = character(),
    lower = numeric(),
    upper = numeric()
  )
  
  for (leaf in 1:length(noderules)) {
    
    rules_tog <- noderules[leaf]
    
    parsed <- strsplit(rules_tog, split = " & ")
    
    leaf_subspace_rules <- tibble(
      leafID = numeric(),
      parameter = character(),
      bound_type = character(),
      bound_value = numeric()
    )
    for (rID in 1:length(parsed[[1]])) {
      
      rule <- parsed[[1]][rID]
      
      if(grepl('<', rule)) {
        
        ruletype_tmp <- 'upper_bd'
        temp_parse <- strsplit(rule, split = ' <= ')
        parameter_tmp <- temp_parse[[1]][1]
        bound_tmp <- as.numeric(temp_parse[[1]][2])
        
        leaf_subspace_rules %<>% bind_rows(tibble(
          leafID = leaf,
          parameter = parameter_tmp,
          bound_type = ruletype_tmp,
          bound_value = bound_tmp
        ))
        
      } else {
        
        ruletype_tmp <- 'lower_bd'
        temp_parse <- strsplit(rule, split = ' >= ')
        parameter_tmp <- temp_parse[[1]][1]
        bound_tmp <- as.numeric(temp_parse[[1]][2])
        
        leaf_subspace_rules %<>% bind_rows(tibble(
          leafID = leaf,
          parameter = parameter_tmp,
          bound_type = ruletype_tmp,
          bound_value = bound_tmp
        ))
        
      }
      
    }
    
    for(param in unique(paramspace$parameter)) {
      
      tmp_fullspace <- paramspace %>%
        filter(parameter == param)
      
      #identify subspace upper bound
      tmp_bounds <- leaf_subspace_rules %>%
        filter(parameter == param, bound_type == 'upper_bd')
      
      # condense to only lowest upper bound
      if (nrow(tmp_bounds) > 0) { # also handles single-bound
        upper_tmp <- min(tmp_bounds$bound_value)
      } else {
        upper_tmp <- tmp_fullspace$maximum[1]
      }
      
      #identify subspace lower bound
      tmp_bounds <- leaf_subspace_rules %>%
        filter(parameter == param, bound_type == 'lower_bd')
      
      # condense to only highest lower bound
      if (nrow(tmp_bounds) > 0) { # also handles single-bound
        lower_tmp <- max(tmp_bounds$bound_value)
      } else {
        lower_tmp <- tmp_fullspace$minimum[1]
      }
      
      leafspace %<>% bind_rows(tibble(
        leafID = leaf,
        parameter = param,
        lower = lower_tmp,
        upper = upper_tmp
      ))
    }
    
  }
  return(leafspace %>%
           mutate(range = upper - lower))
}
