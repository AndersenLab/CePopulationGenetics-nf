#!/usr/bin/env Rscript

require(PopGenome)
require(data.table)
require(tidyr)
require(dplyr)
require(glue)

# args
# 1 - chromosome
# 2 - vcf (${vcf})
# 3 - ancestor (params.anc)
# 4 - window size (params.popgenome_window)
# 5 - window slide (params.popgenome_slide)
# 6 - gff file (params.gff)

# example
# args <- c("I", "some.vcf.gz", "XZ1516", "10000", "1000")

args <- commandArgs(TRUE)

system(glue::glue("echo Initializing PopGenome Parameters"))

chr1 <- c(1,15072434)
chr2 <- c(1,15279421)
chr3 <- c(1,13783801)
chr4 <- c(1,17493829)
chr5 <- c(1,20924180)
chr6 <- c(1,17718942)

chr.lengths <- list(chr1,chr2,chr3,chr4,chr5,chr6)
chroms <- c("I","II","III","IV","V","X")

ANALYSIS_CHROM <- as.character(args[1])
CHROM_START <- chr.lengths[which(chroms == as.character(args[1]))][[1]][1]
CHROM_END <- chr.lengths[which(chroms == as.character(args[1]))][[1]][2]

WINDOW_SIZE <-  as.numeric(args[4])
SLIDE_DISTANCE <- as.numeric(args[5])

OUTGROUP <- as.character(args[3])

system(glue::glue("echo Done Initializing PopGenome Parameters - WindowSize = {WINDOW_SIZE}, StepSize = {SLIDE_DISTANCE}, Whole Population, Chromosome = {ANALYSIS_CHROM}"))

POPGENOME_VCF <- as.character(args[2])
POPGENOME_GFF <- as.character(args[6])

system(glue::glue("echo PopGenome - Reading VCF file")) 
GENOME_OBJECT <- PopGenome::readVCF(
  POPGENOME_VCF, 
  numcols = 10000, 
  tid = ANALYSIS_CHROM, 
  frompos = CHROM_START, 
  topos = CHROM_END, 
  approx = F,
  include.unknown = T)

system(glue::glue("echo PopGenome - Setting Outgroup and Defining Window Size"))

GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)

GENOME_OBJECT <- PopGenome::sliding.window.transform(
  GENOME_OBJECT, 
  width = WINDOW_SIZE, 
  jump = SLIDE_DISTANCE,
  type = 2, 
  whole.data = FALSE
)

system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Detail Stats"))
GENOME_OBJECT <- PopGenome::detail.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Neutrality Stats"))
GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT, detail = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, pi = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Linkage Stats"))
GENOME_OBJECT <- PopGenome::linkage.stats(GENOME_OBJECT, do.ZnS = TRUE, do.WALL = TRUE)

system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

save(GENOME_OBJECT, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Statistics.Rda"))

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Linkage Stats"))

l_df <- data.frame(PopGenome::get.linkage(GENOME_OBJECT)[[1]]) %>%
  tibble::rownames_to_column() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = args[1],
                startWindow = as.numeric(strsplit(rowname,split = " ")[[1]][4]),
                endWindow = as.numeric(strsplit(rowname,split = " ")[[1]][6])) %>%
  dplyr::select(CHROM:endWindow, Wall.B:Kelly.Z_nS) %>%
  tidyr::gather(LinkageStat, value, -(CHROM:endWindow))

save(l_df, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Linkage_Statistics.Rda"))

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Neutrality and Diversity Stats"))

n_df <- data.frame(PopGenome::get.neutrality(GENOME_OBJECT)[[1]]) %>%
  tibble::rownames_to_column() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = args[1],
                startWindow = as.numeric(strsplit(rowname,split = " ")[[1]][4]),
                endWindow = as.numeric(strsplit(rowname,split = " ")[[1]][6])) %>%
  dplyr::select(CHROM:endWindow, Tajima.D:Zeng.E)

d_df <- data.frame(PopGenome::get.diversity(GENOME_OBJECT)[[1]]) %>%
  tibble::rownames_to_column() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = args[1],
                startWindow = as.numeric(strsplit(rowname,split = " ")[[1]][4]),
                endWindow = as.numeric(strsplit(rowname,split = " ")[[1]][6])) %>%
  dplyr::select(CHROM:endWindow, nuc.diversity.within:Pi)

d2_df <- dplyr::left_join(n_df,d_df, by = c("CHROM", "startWindow", "endWindow")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(theta_Watterson = c(GENOME_OBJECT@theta_Watterson)) %>%
  tidyr::gather(statistic, value, -(CHROM:endWindow))

save(d2_df, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Neutrality_Diversity_Statistics.Rda"))