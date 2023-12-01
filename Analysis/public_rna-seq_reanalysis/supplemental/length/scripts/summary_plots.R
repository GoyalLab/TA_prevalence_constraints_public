### Summary plots for length and condition analyses for Mellis et al., 2023
### Created by Madeline E Melzer on 20231120
### Last edit by Madeline E Melzer on 20231128


library(tidyverse)
library(svglite)
library(ggrepel)
library(scales)
#library(ggplot2)
#library(ggbreak)
library(ggbeeswarm)
library(patchwork)
library(tools)

set.seed(23)

plotDirectory = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/transcriptionalcompensation_20231117/plots/"
grn_nitc_path = "/Users/mem3579/Library/CloudStorage/OneDrive-NorthwesternUniversity/Arispe and Goyal Labs/transcriptionalCompensation/grn_nitc/"

#PC (for Madeline's laptop)
plotDirectory = "C:\\Users\\madel\\OneDrive - Northwestern University\\Arispe and Goyal Labs\\transcriptionalCompensation\\transcriptionalcompensation_20231117\\plots\\"
grn_nitc_path = "C:\\Users\\madel\\OneDrive\\Documents\\GitHub\\grn_nitc\\"


allData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path, 'rnaseq/supp_analyses/length/data/allData.csv')), stringsAsFactors = F, header = T))
filteredData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/filteredData.csv')), stringsAsFactors = F, header = T))
mainData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/mainData.csv')), stringsAsFactors = F, header = T))

mechanismComponents <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/nitc_components/data/mechanismComponents_baseMeans.csv')), stringsAsFactors = F, header = T))

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













