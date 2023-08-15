library(dorothea)
library(decoupleR)
library(ggplot2)
library(dplyr)

net <- decoupleR::get_dorothea(levels = c('A', 'B', 'C'))

ko_targets <- c("IRF3", "IRF3", "IRF3", "KLF4", "KLF4", "KLF4", "KLF4", "KLF4", "KLF4", "KLF4", "KLF4", "KLF4", "KLF4", "SOX4", "SOX4", "SOX4", "SOX4", "SOX4", "SOX4", "SOX4", "SOX4", "SOX4", "SOX4", "SP100", "SP100", "SP100", "SP100", "SP100", "SP100", "SP100", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "FOXM1", "STAT3", "E2F1", "E2F1", "E2F1", "E2F1", "FOS", "FOS", "FOS", "SMAD3", "TP53", "TP53", "MYC", "MYC", "STAT1", "STAT1", "STAT1")
paralog_targets <- c("IRF5", "IRF6", "IRF8", "KLF11", "KLF12", "KLF13", "KLF15", "KLF17", "KLF5", "SP5", "SP6", "SP8", "SP9", "CFAP65", "SOX11", "SOX13", "SOX15", "SOX18", "SOX30", "SOX5", "SOX7", "SOX8", "SOX9", "HMG20B", "HMGXB4", "SMARCE1", "SP140", "TFAM", "TOX", "TOX4", "FOXA3", "FOXB1", "FOXD2", "FOXD3", "FOXD4", "FOXD4L1", "FOXD4L5", "FOXD4L6", "FOXE1", "FOXH1", "FOXI3", "FOXJ1", "FOXK2", "FOXL1", "FOXL2", "FOXO6", "FOXR1", "FOXS1", "STAT5A", "E2F2", "E2F4", "E2F5", "E2F8", "BATF", "BATF2", "FOSB", "SMAD6", "TP63", "TP73", "MYCL", "MYCN", "STAT2", "STAT4", "STAT5A")


# create empty dataframe to store results
result_df <- data.frame(
  Value = numeric(),
  stringsAsFactors = FALSE
)

# loop through pairs of symbols
for (i in seq(1, length(ko_targets))) {
  symbol1 <- ko_targets[i]
  symbol2 <- paralog_targets[i]
  
  ko <- filter(net, source == symbol1)
  paralog <- filter(net, source == symbol2)

  res <- inner_join(ko, paralog, by = c('target'))

  result_df <- rbind(result_df, res)
}

# save file to desktop
write.csv(result_df, file = "~/Desktop/paralog_tf_analysis.csv")
