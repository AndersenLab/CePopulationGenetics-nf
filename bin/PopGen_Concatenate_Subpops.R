#!/usr/bin/env Rscript
library(tidyverse)

# PLOT NEUTRALITY STATISTICS
Neutrality_files <- grep(glue::glue("_ND_Statistics"),list.files(), value = T)

neutrality_ls <- list()
for( i in 1:length(Neutrality_files) ){
  
  chr <- strsplit(strsplit(Neutrality_files[i],split = "_")[[1]][1], split="-")[[1]][[2]]
  kpop <- as.numeric(gsub("K", "",strsplit(Neutrality_files[i],split = "_")[[1]][2]))
  if(chr == "23"){
    chr <- "X"
  }
  
  load(Neutrality_files[i])
  neutrality_ls[[i]] <- data.frame(neutrality_df) %>%
    dplyr::mutate(CHROM = chr,
                  K = kpop)
}

neutrality_df <- dplyr::bind_rows(neutrality_ls)

save(neutrality_df, file = glue::glue("Ce_Genome-wide_Neutrality_stats.Rda"))


fst_files <- grep(glue::glue("_FST_Statistics"),list.files(), value = T)

fst_ls <- list()
for( i in 1:length(fst_files)){
  
  chr <- strsplit(strsplit(fst_files[i],split = "_")[[1]][1], split="-")[[1]][[2]]
  kpop <- as.numeric(gsub("K", "",strsplit(Neutrality_files[i],split = "_")[[1]][2]))
  if(chr == "23"){
    chr <- "X"
  }
  load(fst_files[i])
  fst_ls[[i]] <- data.frame(pairFst) %>%
    dplyr::mutate(CHROM = chr,
                  K = kpop)
}

fst_df <- dplyr::bind_rows(fst_ls)

save(fst_df, file = glue::glue("Ce_Genome-wide_Fst_stats.Rda"))
