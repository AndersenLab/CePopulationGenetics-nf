#!/usr/bin/env Rscript
library(tidyverse)

# PLOT NEUTRALITY STATISTICS
Neutrality_files <- grep(glue::glue("Neutrality"),list.files(), value = T)

neutrality_ls <- list()
for( i in 1:length(Neutrality_files) ){
  
  chr <- strsplit(strsplit(Neutrality_files[i],split = "_")[[1]][1], split="-")[[1]][[2]]
  if(chr == "23"){
    chr <- "X"
  }
  
  load(Neutrality_files[i])
  neutrality_ls[[i]] <- data.frame(d2_df) %>%
    dplyr::mutate(CHROM = chr)
}

neutrality_df <- dplyr::bind_rows(neutrality_ls)

save(neutrality_df, file = glue::glue("Ce_Genome-wide_Neutrality_stats.Rda"))

for(neutrality_statistic in unique(neutrality_df$statistic)){
  
  neutrpl <- neutrality_df %>%
    dplyr::filter(statistic == neutrality_statistic) %>%
    ggplot()+
    aes(x = startWindow/1e6, y = value) +
    geom_point(alpha = 0.25 ) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_bw()+
    facet_grid(.~CHROM, scales ="free", space = "free_x") +
    theme(axis.text.x = ggplot2::element_text(size = 16),
          axis.text.y = ggplot2::element_text(size = 16),
          axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
          axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3)) +
    labs(x = "Genomic Position (Mb)",
         y = neutrality_statistic)
  
  ggsave(plot = neutrpl,
         glue::glue("Genome-wide_{neutrality_statistic}.png"),
         width = 12,
         height = 4,
         dpi = 300)
}


Linkage_files <- grep(glue::glue("Linkage"),list.files(), value = T)

linkage_ls <- list()
for( i in 1:length(Linkage_files)){
  
  chr <- strsplit(strsplit(Linkage_files[i],split = "_")[[1]][1], split="-")[[1]][[2]]
  if(chr == "23"){
    chr <- "X"
  }
  load(Linkage_files[i])
  linkage_ls[[i]] <- data.frame(l_df) %>%
    dplyr::mutate(CHROM = chr)
}

linkage_df <- dplyr::bind_rows(linkage_ls)

save(linkage_df, file = glue::glue("Ce_Genome-wide_Linkage_stats.Rda"))

for(linkage_statistic in unique(linkage_df$LinkageStat)){
  
  linkpl <- linkage_df %>%
    dplyr::filter(LinkageStat == linkage_statistic) %>%
    ggplot()+
    aes(x = startWindow/1e6, y = value) +
    geom_point(size = 0.5, alpha = 0.5)+
    theme_bw() +
    facet_grid(.~CHROM, scales ="free", space = "free_x") +
    theme(axis.text.x = ggplot2::element_text(size = 16),
          axis.text.y = ggplot2::element_text(size = 16),
          # legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
          axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3)) +
    labs(x = "Genomic Position (Mb)",
         y = linkage_statistic)
  
  ggsave(plot = linkpl, 
         glue::glue("Genome-wide_{linkage_statistic}.png"),
         width = 12,
         height = 4,
         dpi = 300)
}
