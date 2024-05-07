### Summary plots for length and condition analyses for Mellis et al., 2023
### Created by Madeline E Melzer on 20231120
### Last edit by Madeline E Melzer on 20240503

library(tidyverse)
library(svglite)
library(ggrepel)
library(scales)
library(ggplot2)
#library(ggbreak)
library(ggbeeswarm)
library(patchwork)
library(tools)

set.seed(23)

plotDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/MellisEtAl2023/plots/"
grn_nitc_path = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/MadelineMelzer/MellisEtAl2023/scripts/grn_nitc/"

allData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path, 'rnaseq/supp_analyses/length/data/allData.csv')), stringsAsFactors = F, header = T))
filteredData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/filteredData.csv')), stringsAsFactors = F, header = T))
mainData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/mainData.csv')), stringsAsFactors = F, header = T))

mechanismComponents <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/nitc_components/data/mechanismComponents_baseMeans.csv')), stringsAsFactors = F, header = T))
GTExPairwiseCor <- as_tibble(read.csv(paste0(file.path(grn_nitc_path, 'resub1/Analysis/public_rna-seq_reanalysis/supplemental/GTEx_pariwisecor_FC2.csv')), stringsAsFactors = F, header = T))

############### Condition comparison within one dataset, GSE145653 (-1 and -2) ##############################################################

# Filter and summarize the data
reshaped_data <- allData %>% 
  filter(dataset %in% c('GSE145653-1', 'GSE145653-2')) %>%
  group_by(Gene) %>%
  summarize(`GSE145653-1` = mean(Percent.Upregulated[dataset == 'GSE145653-1']), #these values should all be the same so the average is equivalent
            `GSE145653-2` = mean(Percent.Upregulated[dataset == 'GSE145653-2']),
            pctUpreg_p_value = mean(pctUpreg_p_value[dataset == 'GSE145653-2']) -mean(pctUpreg_p_value[dataset == 'GSE145653-1'])) %>%
  filter(!is.na(`GSE145653-1`) & !is.na(`GSE145653-2`))

# Create the scatter plot
condition_scatter <- ggplot(reshaped_data, aes(x = `GSE145653-1`, y = `GSE145653-2`)) +
  geom_point() +
  geom_text_repel(aes(label = Gene), max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  #scale_color_gradientn(colors = rainbow(3), name = '2-1 %upreg p val')+
  #scale_color_distiller(palette = 'rainbow', name = '%upregulation p.val') +
  xlab('GSE145653-1 fraction upregulated') +
  ylab('GSE145653-2 fraction upregulated') +
  theme_classic()
print(condition_scatter)
#ggsave(condition_scatter, filename = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/condition_scatter.svg", width = 6, height = 4.5)
#ggsave(condition_scatter, filename = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/condition_scatter.png", width = 6, height = 4.5)



### from IAM: count numbers of paralogs per knockout target:
reshaped_data_nPara <- allData %>% 
  dplyr::select(Gene, dataset, Percent.Upregulated) %>%
  filter(dataset %in% c('GSE145653-1', 'GSE145653-2')) %>%
  group_by(Gene, dataset)  %>% 
  mutate(nPara = length(Gene)) %>%
  unique() %>% arrange(Gene) %>%
  pivot_wider(names_from = dataset, values_from = Percent.Upregulated) %>%
  filter(!is.na(`GSE145653-1`),
         !is.na(`GSE145653-2`))

# new version (20240423) of figure 1F that includes number of paralogs per target
condition_scatter_nPara <- ggplot(reshaped_data_nPara, aes(x = `GSE145653-1`, y = `GSE145653-2`)) +
  #geom_point(data = reshaped_data_nPara, aes(x = `GSE145653-1`, y = `GSE145653-2`, size = nPara)) +
  geom_point(aes(color = nPara), size = 2) +  # Use color for nPara
  geom_text_repel(aes(label = Gene), max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_gradient(low = "#ffe4f8", high = "#74003e", name = "Number of Paralogs") +  # Define color gradient
  xlab('GSE145653-1 fraction upregulated') +
  ylab('GSE145653-2 fraction upregulated') +
  guides(color = "none") +
  theme_classic()
print(condition_scatter_nPara)
ggsave(condition_scatter_nPara, filename = paste0(file.path(plotDirectory, 'condition_scatter_nPara.svg')), width = 6, height = 4.5) #20240503MEM
ggsave(condition_scatter_nPara, filename = paste0(file.path(plotDirectory, 'condition_scatter_nPara.png')), width = 6, height = 4.5) #20240503MEM


#### plotting FC2 for each paralog in in condition -1 against condition -2
reshaped_data_conditions <- allData %>% 
  filter(dataset %in% c('GSE145653-1', 'GSE145653-2')) %>%
  select(Gene, Paralog, FC2, dataset) %>%
  group_by(Gene, Paralog) %>%
  pivot_wider(names_from = dataset, values_from = FC2) %>%
  ungroup() %>%
  filter(!is.na(`GSE145653-1`) & !is.na(`GSE145653-2`)) 

pvals <- allData %>%
  filter(dataset == 'GSE145653-1') %>%
  select(Gene, Paralog, padj) %>%
  mutate(padj_transform = -log10(padj))

reshaped_data_conditions <- left_join(reshaped_data_conditions, pvals, by = c("Gene", "Paralog"))

max_padj_transform <- max(reshaped_data_conditions$padj_transform, na.rm = TRUE)

reshaped_data_conditions_noNA = reshaped_data_conditions %>% filter(is.na(padj) == FALSE)

condition_scatter_FC2 <- ggplot(reshaped_data_conditions_noNA, aes(x = `GSE145653-1`, y = `GSE145653-2`, color = padj_transform)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_gradientn(
    colors = c("#e3cfff", "#0000ff", "#0b003f", "#0b003f"),  # Ensuring dark color for high values
    values = rescale(c(0, 14.99, 15, max_padj_transform)),  # Ensuring transition ends at 15
    name = "-log10(padj)",
    limits = c(0, max_padj_transform),  # Using full range of data for plotting
    breaks = c(0, 5, 10, 15),  # Define breaks to show in the legend
    labels = c("0", "5", "10", "15+"),  # Define labels corresponding to breaks
    guide = guide_colorbar(title = "-log10(padj)", barwidth = 2, barheight = 6, ticks = TRUE)  # Color bar settings
  ) +
  xlab('GSE145653-1 FC2') +
  ylab('GSE145653-2 FC2') +
  theme_classic()
print(condition_scatter_FC2)
#ggsave(condition_scatter_FC2, filename = paste0(file.path(plotDirectory, 'condition_scatter_FC2.svg')), width = 6, height = 4.5) #20240502MEM
#ggsave(condition_scatter_FC2, filename = paste0(file.path(plotDirectory, 'condition_scatter_FC2.png')), width = 6, height = 4.5) #20240502MEM


# one plot per gene

reshaped_data_conditions_hits = reshaped_data_conditions_noNA %>% filter(Gene %in% c("Lmna", "Macf1", "Myc", "Nes", "Nsd1", "Trim71", "Zfp281"))

condition_scatter_FC2 <- ggplot(reshaped_data_conditions_noNA, aes(x = `GSE145653-1`, y = `GSE145653-2`, color = padj_transform)) +
  geom_point() +
  #geom_text_repel(aes(label = Paralog), max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ Gene, scales = "free") +  # Create a separate plot for each gene
  scale_color_gradientn(
    colors = c("#e3cfff", "#0000ff", "#0b003f", "#0b003f"),  # Ensuring dark color for high values
    values = rescale(c(0, 14.99, 15, max_padj_transform)),  # Ensuring transition ends at 15
    name = "-log10(padj)",
    limits = c(0, max_padj_transform),  # Using full range of data for plotting
    breaks = c(0, 5, 10, 15),  # Define breaks to show in the legend
    labels = c("0", "5", "10", "15+"),  # Define labels corresponding to breaks
    guide = FALSE
  ) +
  xlab('GSE145653-1 FC2') +
  ylab('GSE145653-2 FC2') +
  theme_classic()
print(condition_scatter_FC2)
#ggsave(condition_scatter_FC2, filename = paste0(file.path(plotDirectory, 'condition_scatter_FC2_all.svg')), width = 7.5, height = 8) #20240502MEM, saved for hits (6 x 4.5) and all individuals
#ggsave(condition_scatter_FC2, filename = paste0(file.path(plotDirectory, 'condition_scatter_FC2_all.png')), width = 7.5, height = 8) #20240502MEM, saved for hits (6 x 4.5) and all individuals


##stats, 20240502
cor(reshaped_data_conditions_noNA$"GSE145653-1", reshaped_data_conditions_noNA$"GSE145653-2", method = "spearman") # 0.6772401, 20240502


############### gene length analysis #########################################################################################################

##### Coloring scatter plot by NITC or non-NITC gene (this is what I sent to Yogesh to make final figure)

# Define the list of specific genes
specific_genes <- c('ASH1L', 'Actb', 'Apc', 'KLF6', 'Lmna', 'Macf1', 'Myc', 'Nes', 'Nsd1', 'RB1', 'RUNX3', 'SP1', 'TLR4', 'Trim71', 'ZEB2', 'Zfp281')

# Reshape data
reshaped_data <- mainData %>%
  group_by(Gene, dataset) %>%
  filter(!is.na(ko_length)) %>%
  summarize(ko_length = mean(ko_length, na.rm = TRUE),
            Percent.Upregulated = mean(Percent.Upregulated, na.rm = TRUE),
            pctUpreg_p_value = mean(pctUpreg_p_value)) %>%
  mutate(gene_category = if_else(Gene %in% specific_genes, "NITC", "non-NITC"))

# Create the scatter plot with the new color aesthetic
length_scatter <- ggplot(reshaped_data, aes(x = ko_length, y = Percent.Upregulated, color = gene_category)) +
  geom_point() +
  geom_text_repel(aes(label = Gene), max.overlaps = 8) +
  #geom_abline(slope = 1, intercept = 0, linetype = 2) +
  #geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("non-NITC" = "grey", "NITC" = "black"), name = 'Gene Category') +
  xlab('Length of CRISPR KO Gene') +
  ylab('% Upregulation') +
  #xlim(0, 600000) +
  #annotate("text", x = Inf, y = -Inf, label = "Linear Trendline: y ~ x", hjust = 1.5, vjust = -11, size = 4, color = "blue") +
  theme_classic()
print(length_scatter)
ggsave(length_scatter, filename = "C:/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/length_scatter_cut.svg", width = 6, height = 4.5)
ggsave(length_scatter, filename = "C:/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/length_scatter_cut.png", width = 6, height = 4.5)



##### Plotting violin plot of length vs NITC or non-NITC genes

reshaped_data <- reshaped_data %>% group_by(Gene) %>% summarize(ko_length = mean(ko_length)) %>%
  mutate(gene_category = if_else(Gene %in% specific_genes, "NITC", "non-NITC"))

length_violin = ggplot(reshaped_data, aes(x = gene_category, y = ko_length)) +
  geom_violin(aes(fill = gene_category), trim = FALSE, color = NA) + 
  geom_point(color = "black", position = position_jitter(width = 0.2)) +  
  xlab("Gene Category") +
  ylab("KO Length") +
  ggtitle("KO Length by Gene Category") +
  theme_classic()
length_violin
#ggsave(length_violin, filename = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/length_violin_NITC_vs_non.svg", width = 6, height = 4.5)
#ggsave(length_violin, filename = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/length_violin_NITC_vs_non.png", width = 6, height = 4.5)








############### NITC mechanistically-relevant genes: components of NDM-regulation and COMPASS complex #########################################################################################################

# Define the list of specific genes
specific_genes <- c('ASH1L', 'Actb', 'Apc', 'KLF6', 'Lmna', 'Macf1', 'Myc', 'Nes', 'Nsd1', 'RB1', 'RUNX3', 'SP1', 'TLR4', 'Trim71', 'ZEB2', 'Zfp281') 

crispr_targets <- c('ASH1L', 'Actb', 'Actg1', 'Aff3', 'Apc', 'BUB1B', 'CMTR1', 'Cdk8', 'Csnk1a1',
                      'Etv5', 'FOSL1', 'Fbxw7', 'Fermt2', 'Furin', 'Hprt', 'IKZF5', 'INO80', 'JUNB',
                      'Jarid2', 'Jmjd1c', 'KDM1A', 'KLF6', 'KMT2A', 'L3mbtl3', 'LCK', 'Lmna', 'MBD2',
                      'MTF1', 'Macf1', 'Mbd3', 'Msi2', 'Myc', 'Myo10', 'NSD1', 'Nes', 'Nmt1', 'Nsd1',
                      'PLK1', 'POLK', 'PRPF4B', 'Pten', 'Pum1', 'RB1', 'REL', 'RUNX3', 'Raf1', 'RhoC',
                      'SIX6', 'SMARCA4', 'SMC3', 'SOX10', 'SP1', 'SRF', 'STAG2', 'STK11', 'Smg7', 'TFAM',
                      'TLR4', 'TP53', 'Tcf7l1', 'Tet1', 'Trim71', 'UBR5', 'Usp7', 'Usp9x', 'VEZF1',
                      'WDR7', 'ZEB2', 'Zfp281', 'Zfp423', 'p50', 'p52', 'p65')


# keeping mouse and human genes the same format so they are plotted together (ONLY RUN IF YOU WANT MOUSE-HUMAN MIX)
mechanismComponents <- mechanismComponents %>%
  mutate(component_gene = tolower(as.character(component_gene))) %>%
  mutate(component_gene = paste0(toupper(substring(component_gene, 1, 1)), 
                                 substring(component_gene, 2)))


# Reshape data
reshaped_data <- mechanismComponents %>%
  filter(dataset == 'GSE145653-1') %>% # only plotting one dataset at a time
  #filter(!is.na(Percent.Upregulated)) %>%
  mutate(gene_category = if_else(ko_gene %in% specific_genes, "NITC", "non-NITC")) %>%
  filter(ko_gene %in% crispr_targets)

component_plots = list()


component_genes = unique(mechanismComponents$component_gene) #only want to make one plot for each mechanism component gene

for (gene in component_genes) {
  component_gene_data <- reshaped_data %>% 
    filter(component_gene == gene) %>%
    mutate(gene_category = factor(gene_category, levels = c("non-NITC", "NITC")))
  
  p <- ggplot(component_gene_data, aes(x = gene_category, y = baseMean, fill = gene_category, color = gene_category)) +
    geom_beeswarm(aes(color = gene_category), size = 1, cex = 1.5) +
    #geom_violin(trim = FALSE) +
    ggtitle(paste("Gene:", gene)) + #for the all-human dataset GSE151825, should be capitalized, so case matches species convention
    xlab('Gene Category') +
    ylab('Base Mean') +
    scale_color_manual(values = c("non-NITC" = "grey", "NITC" = "black")) +
    scale_fill_manual(values = c("non-NITC" = "grey", "NITC" = "black")) +
    theme(legend.position = "none") +
    theme_classic()
  
  # Store the plot in the list
  component_plots[[gene]] <- p
}

plot_grid <- reduce(component_plots[13:24], `+`) + 
  plot_layout(ncol = 4) &
  theme(legend.position = "none")

print(plot_grid)
#ggsave(plot_grid, filename = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/nitcComponents_beeswarm_GSE145653-1.svg", width = 6, height = 4.5)
#ggsave(plot_grid, filename = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/nitcComponents_violin.png", width = 6, height = 4.5)

# Extract unique gene names for NITC and non-NITC
unique_genes_NITC <- reshaped_data %>%
  filter(gene_category == "NITC") %>%
  distinct(ko_gene)

unique_genes_non_NITC <- reshaped_data %>%
  filter(gene_category == "non-NITC") %>%
  distinct(ko_gene)

# Print the unique gene names
print("Unique genes in NITC category:")
print(unique_genes_NITC)

print("Unique genes in non-NITC category:")
print(unique_genes_non_NITC)



############### whether coexpression is predictive of paralog log2foldchange after target knockout #########################################################################################################
# from IAM analysis, 20240501

GTExPairwiseCor = GTExPairwiseCor %>%
  mutate(padj_transform = -log10(padj))

max_padj_transform <- max(GTExPairwiseCor$padj_transform, na.rm = TRUE)

GTExPairwiseCor_noNA = GTExPairwiseCor %>% filter(is.na(padj) == FALSE)

GTExPairwiseCor_plot <- ggplot(GTExPairwiseCor_noNA, aes(x = cor_spearman, y = FC2, color = padj_transform)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_gradientn(
    colors = c("#e3cfff", "#0000ff", "#0b003f", "#0b003f"),  # Ensuring dark color for high values
    values = rescale(c(0, 9.99, 10, max_padj_transform)),  # Ensuring transition ends at 15
    name = "-log10(padj)",
    limits = c(0, max_padj_transform),  # Using full range of data for plotting
    breaks = c(0, 5, 10),  # Define breaks to show in the legend
    labels = c("0", "5", "10+"),  # Define labels corresponding to breaks
    #guide = FALSE
    guide = guide_colorbar(title = "-log10(padj)", barwidth = 2, barheight = 6, ticks = TRUE)  # Color bar settings
  ) +
  xlab("Spearman's correlation") +
  ylab("FC2") +
  theme_classic()
print(GTExPairwiseCor_plot)
ggsave(GTExPairwiseCor_plot, filename = paste0(file.path(plotDirectory, 'GTExPairwiseCor.svg')), width = 6, height = 4.5) #20240502MEM
ggsave(GTExPairwiseCor_plot, filename = paste0(file.path(plotDirectory, 'GTExPairwiseCor.png')), width = 6, height = 4.5) #20240502MEM


# for scalebar ONLY
GTExPairwiseCor_plot <- ggplot(GTExPairwiseCor_noNA, aes(x = cor_spearman, y = FC2, color = padj_transform)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_gradientn(
    colors = c("#e3cfff", "#0000ff"),  # Ensuring dark color for high values
    values = rescale(c(0, 10)),  # Ensuring transition ends at 15
    name = "-log10(padj)",
    limits = c(0, 10),  # Using full range of data for plotting
    breaks = c(0, 5, 10),  # Define breaks to show in the legend
    labels = c("0", "5", "10+"),  # Define labels corresponding to breaks
    guide = guide_colorbar(title = "-log10(padj)", barwidth = 2, barheight = 6, ticks = TRUE)  # Color bar settings
  ) +
  xlab("Spearman's correlation") +
  ylab("FC2") +
  theme_classic()
print(GTExPairwiseCor_plot)
#ggsave(GTExPairwiseCor_plot, filename = paste0(file.path(plotDirectory, 'GTExPairwiseCor_scalebar.svg')), width = 6, height = 4.5) #20240502MEM


