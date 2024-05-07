#TADKB

library(tidyverse)
library(magrittr)
library(biomaRt)
library(gprofiler2)
# Load TAD coordinates - load all DI files at 50kb resolution (per webtool); use any method, union
tadkb_files_50kb_DI_anymethod  <- tibble(filen = list.files('~/Documents/TAD_annotations/TADs/')) %>%
  filter(grepl('50kb', filen), grepl('DI', filen)) %>%
  mutate(species = ifelse(grepl('_CH12', filen) | 
                            grepl('_ES_', filen) | 
                            grepl('_NPC_', filen) | 
                            grepl('_CN_', filen), 'mouse', 'human')) %>%
  filter(species == 'human') %>%
  dplyr::select(-species)

all_tads <- tibble(
  chr = character(),
  start = numeric(),
  end = numeric(),
  TAD_fn = character()
)
for(tf in 1:nrow(tadkb_files_50kb_DI_anymethod)) {
  tmp_tads <- as_tibble(read.table(paste0('~/Documents/TAD_annotations/TADs/', tadkb_files_50kb_DI_anymethod$filen[tf]))) %>%
    mutate(TAD_fn = tadkb_files_50kb_DI_anymethod$filen[tf]) %>%
    dplyr::rename(chr = V1,
                  start = V2,
                  end = V3)
  all_tads %<>% bind_rows(tmp_tads)
}

# Get gene loci from biomaRt (TSS)
mainData <- as_tibble(read.csv(paste0(file.path(grn_nitc_path,'rnaseq/supp_analyses/length/data/mainData.csv')), stringsAsFactors = F, header = T))

mainData_human <- mainData %>% filter(species == 'human')
conv_HGNCs_targs <- gconvert(
  query = mainData_human$Gene,
  organism = "hsapiens",
  target = "HGNC",
  numeric_ns = "",
  mthreshold = 1,
  filter_na = FALSE
)
conv_HGNCs_paras <- gconvert(
  query = mainData_human$Paralog,
  organism = "hsapiens",
  target = "HGNC",
  numeric_ns = "",
  mthreshold = 1,
  filter_na = FALSE
)

mainData_human$target_HGNC <- conv_HGNCs_targs$target
mainData_human$paralog_HGNC <- conv_HGNCs_paras$target

mainData_filt <- mainData_human %>% filter(!is.na(target_HGNC), !is.na(paralog_HGNC))

targets <- unique(mainData_human$target_HGNC)
paralogs <- unique(mainData_human$paralog_HGNC)

mart <- useMart("ensembl", "hsapiens_gene_ensembl", 'https://grch37.ensembl.org')
att <- listAttributes(mart)
grep("transcript", att$name, value=TRUE)

target_TSS <- as_tibble(getBM(attributes=c("transcription_start_site", "chromosome_name",
                                 "transcript_start", "transcript_end",
                                 "strand",  "ensembl_gene_id",
                                 "ensembl_transcript_id", "external_gene_name","hgnc_symbol"),
                    filters="hgnc_symbol", 
                    values=targets, 
                    mart=mart)) %>%
  dplyr::rename(target_HGNC = hgnc_symbol) %>%
  mutate(chr = paste0('chr', chromosome_name))
paralog_TSS <-  as_tibble(getBM(attributes=c("transcription_start_site", "chromosome_name",
                                             "transcript_start", "transcript_end",
                                             "strand",  "ensembl_gene_id",
                                             "ensembl_transcript_id", "external_gene_name","hgnc_symbol"),
                                filters="hgnc_symbol", 
                                values=paralogs, 
                                mart=mart)) %>%
  dplyr::rename(paralog_HGNC = hgnc_symbol) %>%
  mutate(chr = paste0('chr', chromosome_name))

mainData_pared <- mainData_human %>%
  mutate(isUpreg = FC2>0.5 & padj<0.05) %>%
  dplyr::select(Gene, target_HGNC, Paralog, paralog_HGNC, isUpreg)
# Check overlap with TAD coordinate
mainData_pared_out <- tibble(
  Gene = character(),
  target_HGNC = character(),
  Paralog = character(),
  paralog_HGNC = character(),
  isUpreg = logical(),
  inAnyTad = logical()
)
for (pr in 1:nrow(mainData_human)) {
  
  tmp_mainData <- mainData_pared[pr,]
  
  mainDat_tmp_out <- tmp_mainData 
  
  targ_HGNC <- unique(tmp_mainData$target_HGNC)
  
  tmp_targ_TSS <- target_TSS %>% dplyr::select("target_HGNC", "chr", "transcription_start_site") %>%
    filter(target_HGNC == unique(tmp_mainData$target_HGNC)) %>%
    dplyr::rename(chr_target = chr,
                  transcript_start_target = transcription_start_site)
  tmp_para_TSS <- paralog_TSS %>% dplyr::select("paralog_HGNC", "chr", "transcription_start_site") %>%
    filter(paralog_HGNC == unique(tmp_mainData$paralog_HGNC)) %>%
    dplyr::rename(chr_paralog = chr,
                  transcript_start_paralog = transcription_start_site) %>%
    filter(chr_paralog %in% unique(tmp_targ_TSS$chr_target))
  
  if(nrow(tmp_para_TSS) == 0){
    mainDat_tmp_out$inAnyTad = FALSE
    mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
  } else {
    
    tmp_TADs <- all_tads %>% 
      filter(chr %in% unique(tmp_targ_TSS$chr_target)) 
    tmp_TADs2 <- tibble(
      chr = character(),
      start = numeric(),
      end = numeric(),
      TAD_fn = character(),
      contains_targ = logical()
    )
    for (td in 1:nrow(tmp_TADs)) {
      
      tmp_TADs_row <- tmp_TADs[td,]
      targTADs_tmp <- tmp_TADs_row %>%
        mutate(contains_targ = sum(tmp_targ_TSS$transcript_start_target > start & tmp_targ_TSS$transcript_start_target < end)>0)
      tmp_TADs2 %<>% bind_rows(targTADs_tmp)
      
    }
    
    tmp_TADs_withtarg <- tmp_TADs2 %>%
      filter(contains_targ == TRUE)
    
    if(nrow(tmp_TADs_withtarg) == 0){
      mainDat_tmp_out$inAnyTad = FALSE
      mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
    } else {
      
      min_TAD_start = min(tmp_TADs_withtarg$start)
      max_TAD_end = max(tmp_TADs_withtarg$end)
      
      tmp_para_TSS_inTADrange <- tmp_para_TSS %>%
        filter(transcript_start_paralog > min_TAD_start, transcript_start_paralog < max_TAD_end)
      
      if(nrow(tmp_para_TSS_inTADrange) == 0){
        mainDat_tmp_out$inAnyTad = FALSE
        mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
      } else {
        
        paras_in_tads <- tibble(
          paralog_HGNC = character(),
          inTad = logical(),
          TR = numeric()
        )
        for (tr in 1:nrow(tmp_TADs_withtarg)) {
          
          tad_tmp <- tmp_TADs_withtarg[tr,]
          
          paras_in_tads %<>% bind_rows(tmp_para_TSS_inTADrange %>%
                                         mutate(inTad_tss = transcript_start_paralog > tad_tmp$start & transcript_start_paralog < tad_tmp$end) %>%
                                         group_by(paralog_HGNC) %>%
                                         summarise(inTad = sum(inTad_tss) > 0) %>%
                                         mutate(TR = tr))
          
        }
        paras_in_anytad <- paras_in_tads %>%
          group_by(paralog_HGNC) %>%
          summarise(inAnyTad = sum(inTad) > 0)
        
        mainDat_tmp_out %<>% left_join(paras_in_anytad, by = 'paralog_HGNC')
        mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
      }
    }
  }
}
fisher.test(mainData_pared_out$isUpreg, mainData_pared_out$inAnyTad)

# consider all genes, not just paralogs. Lots more data, so will need to iterate through targets and load respective datasets

# dataset_meta

unique_sets <- unique(mainData_human$GEO_ID)

# targ_anygene_sharestad_all <- list()
# tad_sub_targ_total_anygene <- list()
# anygenes_tad_summaryshare <- list()
TAD_FET_anygene_pvals <- tibble(sampleKO = character(),
       target_HGNC = character(),
       GEO_ID = character(),
       FET_pval = numeric())
for (ds in unique_sets) {
  
  cat('Loading', ds, 'and converting gene IDs...')
  allDES_tmp <- read.csv(paste0(grn_nitc_path, '/rnaseq/deseq_files/', ds, '/differentialExpression_DESeq_allTargets.csv'), stringsAsFactors = F, header = T)
  
  tmp_main <- mainData %>%
    filter(GEO_ID == ds,
           !is.na(target_HGNC))
  
  unique_targs_in_ds <- unique(tmp_main$target_HGNC)
  
  target_TSS <- as_tibble(getBM(attributes=c("transcription_start_site", "chromosome_name",
                                             "transcript_start", "transcript_end",
                                             "strand",  "ensembl_gene_id",
                                             "ensembl_transcript_id", "external_gene_name","hgnc_symbol"),
                                filters="hgnc_symbol", 
                                values=unique_targs_in_ds, 
                                mart=mart)) %>%
    dplyr::rename(target_HGNC = hgnc_symbol) %>%
    mutate(chr = paste0('chr', chromosome_name))
  
  # if (ds != 'GSE151825') {
  conv_HGNCs_anygene_tmp <- gconvert(
    query = unique(allDES_tmp$gene_name),
    organism = "hsapiens",
    target = "HGNC",
    numeric_ns = "",
    mthreshold = 1,
    filter_na = FALSE
  )
  
  allDES_tmp %<>%
    left_join(conv_HGNCs_anygene_tmp %>% 
                dplyr::select(input, target) %>%
                dplyr::rename(gene_name = input,
                              othergene_HGNC = target), by = 'gene_name') %>%
    filter(!is.na(othergene_HGNC))
  
  # allDES_tmp$othergene_HGNC <- conv_HGNCs_anygene_tmp$target
  
  othergene_TSS <-  as_tibble(getBM(attributes=c("transcription_start_site", "chromosome_name","ensembl_gene_id",
                                                 "ensembl_transcript_id", "external_gene_name","hgnc_symbol"),
                                    filters="hgnc_symbol", 
                                    values=unique(allDES_tmp$othergene_HGNC), 
                                    mart=mart)) %>%
    dplyr::rename(othergene_HGNC = hgnc_symbol) %>%
    mutate(chr = paste0('chr', chromosome_name))
  
  # }
  
  for (gs in unique_targs_in_ds) {
    
    if (ds == 'GSE151825') {
      cat('Working on', gs, '\n')
      
      allDES_tmp1 <- allDES_tmp %>%
        filter(sampleKO == gs)

    }
    
    conv_tmp <- conv_HGNCs_targs %>%
      filter(target == gs) %>%
      dplyr::select(input, target) %>% unique()
    
    gs_hgnc <- gs
    
    if (ds == "GSE151825") {
      mainData_tmp <- allDES_tmp1 %>%
        filter(sampleKO == conv_tmp$input) %>%
        mutate(target_HGNC = gs)
    } else {
      mainData_tmp <- allDES_tmp %>%
        filter(sampleKO == conv_tmp$input) %>%
        mutate(target_HGNC = gs)
    }
    
    mainData_tmp_pared <- as_tibble(mainData_tmp) %>%
      mutate(isUpreg = log2FoldChange>0.5 & padj<0.05) %>%
      dplyr::select(sampleKO, target_HGNC, gene_name, othergene_HGNC, isUpreg)
    
    mainDat_tmp_out <- mainData_tmp_pared %>%
      mutate(GEO_ID = ds)
    
    # targ_HGNC <- unique(tmp_mainData$target_HGNC)
    
    tmp_targ_TSS <- target_TSS %>% dplyr::select("target_HGNC", "chr", "transcription_start_site") %>%
      filter(target_HGNC == unique(mainData_tmp_pared$target_HGNC)) %>%
      dplyr::rename(chr_target = chr,
                    transcript_start_target = transcription_start_site)
    
    tmp_other_TSS <- othergene_TSS %>% dplyr::select("othergene_HGNC", "chr", "transcription_start_site") %>%
      # filter(othergene_HGNC == unique(tmp_mainData$paralog_HGNC)) %>%
      dplyr::rename(chr_paralog = chr,
                    transcript_start_othergene = transcription_start_site) %>%
      filter(chr_paralog %in% unique(tmp_targ_TSS$chr_target))
    
    if(nrow(tmp_other_TSS) == 0){
      mainDat_tmp_out$inAnyTad = FALSE
      mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
    } else {
      
      tmp_TADs <- all_tads %>% 
        filter(chr %in% unique(tmp_targ_TSS$chr_target)) 
      tmp_TADs2 <- tibble(
        chr = character(),
        start = numeric(),
        end = numeric(),
        TAD_fn = character(),
        contains_targ = logical()
      )
      for (td in 1:nrow(tmp_TADs)) {
        
        tmp_TADs_row <- tmp_TADs[td,]
        targTADs_tmp <- tmp_TADs_row %>%
          mutate(contains_targ = sum(tmp_targ_TSS$transcript_start_target > start & tmp_targ_TSS$transcript_start_target < end)>0)
        tmp_TADs2 %<>% bind_rows(targTADs_tmp)
        
      }
      
      tmp_TADs_withtarg <- tmp_TADs2 %>%
        filter(contains_targ == TRUE)
      
      if(nrow(tmp_TADs_withtarg) == 0){
        mainDat_tmp_out$inAnyTad = FALSE
        mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
      } else {
        
        min_TAD_start = min(tmp_TADs_withtarg$start)
        max_TAD_end = max(tmp_TADs_withtarg$end)
        
        tmp_others_TSS_inTADrange <- tmp_other_TSS %>%
          filter(transcript_start_othergene > min_TAD_start, transcript_start_othergene < max_TAD_end)
        
        if(nrow(tmp_others_TSS_inTADrange) == 0){
          mainDat_tmp_out$inAnyTad = FALSE
          mainData_pared_out %<>% bind_rows(mainDat_tmp_out)
        } else {
          
          others_in_tads <- tibble(
            othergene_HGNC = character(),
            inTad = logical(),
            TR = numeric()
          )
          for (tr in 1:nrow(tmp_TADs_withtarg)) {
            
            tad_tmp <- tmp_TADs_withtarg[tr,]
            
            others_in_tads %<>% bind_rows(tmp_others_TSS_inTADrange %>%
                                           mutate(inTad_tss = transcript_start_othergene > tad_tmp$start & transcript_start_othergene < tad_tmp$end) %>%
                                           group_by(othergene_HGNC) %>%
                                           summarise(inTad = sum(inTad_tss) > 0) %>%
                                           mutate(TR = tr))
            
          }
          others_in_anytad <- others_in_tads %>%
            group_by(othergene_HGNC) %>%
            summarise(inAnyTad = sum(inTad) > 0)
          
          mainDat_tmp_out %<>% left_join(others_in_anytad, by = 'othergene_HGNC')
          mainDat_tmp_out$inAnyTad %<>% replace_na(FALSE)
          
          
          
        }
      }
    }
    if(length(unique(mainDat_tmp_out$inAnyTad))>1){
      tmp_FET <- fisher.test(mainDat_tmp_out$isUpreg, mainDat_tmp_out$inAnyTad)
      
      TAD_FET_anygene_pvals %<>% bind_rows(tibble(sampleKO = conv_tmp$input[1],
                                                  target_HGNC = gs,
                                                  GEO_ID = ds,
                                                  FET_pval = tmp_FET$p.value))
    } else {
      TAD_FET_anygene_pvals %<>% bind_rows(tibble(sampleKO = conv_tmp$input[1],
                                                  target_HGNC = gs,
                                                  GEO_ID = ds,
                                                  FET_pval = NA))
    }
  }
}
TAD_FET_anygene_pvals$p.adj <- p.adjust(TAD_FET_anygene_pvals$FET_pval, method = 'BH')
write.csv(TAD_FET_anygene_pvals, file = paste0(grn_nitc_path, 'rnaseq/de_analysis/summary_plots_drafts/tadkb_anygene_FET.csv'), quote = F, row.names = F)
