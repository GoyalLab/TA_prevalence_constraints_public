library(dorothea)
library(decoupleR)
library(ggplot2)
library(dplyr)

# download data
human <- decoupleR::get_dorothea(levels = c("A", "B", "C"), organism = "human")
mouse <- decoupleR::get_dorothea(levels = c("A", "B", "C"), organism = "mouse")
net <- rbind(human, mouse)

# Hits --------------
ko_targets <- c("KLF6", "RUNX3", "SP1", "ZEB2", "Myc", "Myc", "Zfp281", "Zfp281", "Zfp281", "Zfp281", "Zfp281", "Zfp281", "Zfp281", "Zfp281", "Zfp281")
paralog_targets <- c("KLF7", "RUNX1", "SP4", "ZEB1", "Mycl", "Mycn", "Zfp276", "Zbtb16", "Zfp692", "Zbtb45", "Snai1", "Prdm6", "Zfp296", "Zfp653", "Zbtb32")

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

  res <- inner_join(ko, paralog, by = c("target"))

  result_df <- rbind(result_df, res)
}

# save file to desktop
write.csv(result_df, file = "~/Desktop/hits_tf_analysis.csv")


# Non-Hits --------------
ko_targets <- c("FOSL1", "FOSL1", "FOSL1", "FOSL1", "FOSL1", "FOSL1", "JUNB", "JUNB", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "MTF1", "SOX10", "SOX10", "SOX10", "SOX10", "SOX10", "SOX10", "SOX10", "SRF", "SRF", "SRF", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "VEZF1", "Aff3", "Aff3", "Aff3", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "Etv5", "L3mbtl3", "Mbd3", "Mbd3", "Tcf7l1", "Tcf7l1", "Tcf7l1", "Zfp423", "Zfp423", "Zfp423", "Zfp423", "Zfp423", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "IKZF5", "REL", "REL", "REL", "REL", "SIX6", "SIX6", "SIX6", "SIX6", "SIX6", "TP53", "TP53", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "TFAM", "MBD2", "MBD2")
paralog_targets <- c("JDP2", "FOS", "ATF3", "FOSB", "BATF", "FOSL2", "JUND", "JUN", "VEZF1", "ZNF384", "ZNF76", "ZBTB20", "MAZ", "PATZ1", "ZFP91", "SNAI1", "ZBTB1", "PRDM10", "HIC2", "ZBTB45", "ZNF276", "ZBTB42", "ZNF148", "SNAI2", "HIC1", "ZNF575", "ZNF281", "ZNF653", "ZNF143", "ZBTB18", "ZNF362", "ZNF692", "SOX8", "SOX5", "SOX13", "SOX12", "SOX6", "SOX4", "SOX15", "MEF2D", "MEF2C", "MEF2A", "ZNF76", "ZBTB20", "PATZ1", "ZNF384", "MAZ", "PRDM10", "ZFP91", "ZBTB42", "ZBTB1", "HIC2", "SNAI1", "ZNF148", "ZBTB45", "ZNF281", "HIC1", "ZNF276", "ZNF575", "ZNF653", "ZNF143", "ZBTB18", "ZNF362", "ZNF692", "MTF1", "SNAI2", "Aff1", "Aff2", "Aff4", "Elf1", "Elk3", "Spib", "Gabpa", "Spic", "Ets1", "Etv3", "Etv1", "Elk1", "Erf", "Ets2", "Elf3", "Etv4", "Elk4", "Elf4", "Etv2", "Etv6", "Spi1", "Elf2", "L3mbtl1", "Mbd1", "Mbd2", "Tcf7", "Tcf7l2", "Lef1", "Zfp597", "Zfp521", "Zfp445", "Zfp786", "Zbtb39", "ZNF587B", "ZNF586", "ZNF837", "ZNF570", "ZNF549", "E4F1", "ZNF773", "ZNF304", "ZNF584", "ZNF583", "ZNF79", "ZNF562", "ZNF561", "ZSCAN2", "ZNF691", "ZNF548", "ZNF480", "ZSCAN20", "ZNF134", "ZNF610", "ZNF552", "ZNF660", "ZNF793", "ZFAT", "ZNF792", "ZNF419", "ZNF711", "ZFX", "ZNF8", "PRDM15", "ZBTB11", "ZNF551", "ZNF772", "ZFY", "ZNF154", "RELB", "NFKB1", "NFKB2", "RELA", "SIX2", "SIX3", "SIX4", "SIX1", "SIX5", "TP73", "TP63", "HMG20B", "HMGB3", "HMGB2", "TOX", "TOX2", "SP100", "TOX3", "SP140L", "HMGB1", "HMGXB4", "UBTF", "SMARCE1", "SP110", "SSRP1", "HMG20A", "TOX4", "MBD1", "MBD3")

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

  res <- inner_join(ko, paralog, by = c("target"))

  result_df <- rbind(result_df, res)
}

# save file to desktop
write.csv(result_df, file = "~/Desktop/non-hits_tf_analysis.csv")
