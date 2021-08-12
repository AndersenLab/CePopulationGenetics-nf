#!/usr/bin/env Rscript

# arguments
# 1 - Chromosome
# 2 - Ancestor: ${params.anc}
# 3 - VCF: ${vcf}
# 4 - GFF: ${params.gff}
# 5 - Strain population text file
# 6 - K size
# 7 - window size
# 8 - slide distance

# testing
# args <- c("V", "XZ1516", "Ce330_annotated.vcf.gz", "/projects/b1059/projects/Stefan/CePopulationGenetics-nf/input_files/Popgenome_files/WS245_exons.gff", "LD_0.8_MAF_0.025.10.Q",  "10", "10000", "1000")

args <- commandArgs(trailingOnly = TRUE)

require(PopGenome)
require(data.table)
require(tidyverse)

system(glue::glue("echo Initializing PopGenome Parameters"))

chr1 <- c(1,15072434)
chr2 <- c(1,15279421)
chr3 <- c(1,13783801)
chr4 <- c(1,17493829)
chr5 <- c(1,20924180)
chr6 <- c(1,17718942)

chr.lengths <- list(chr1,chr2,chr3,chr4,chr5,chr6)
chroms <- c("I","II","III","IV","V","X")

system(glue::glue("echo Done Initializing PopGenome Parameters "))

ANALYSIS_CHROM <- args[1]
CHROM_START <- chr.lengths[which(chroms == ANALYSIS_CHROM)][[1]][1]
CHROM_END <- chr.lengths[which(chroms == ANALYSIS_CHROM)][[1]][2]

OUTGROUP <- args[2]
POPGENOME_VCF <- args[3]
POPGENOME_GFF <- args[4]

WINDOW_SIZE <-  as.numeric(args[7])
SLIDE_DISTANCE <- as.numeric(args[8])

####################### - set up populations
system(glue::glue("bcftools query -l {POPGENOME_VCF} | sort > sample_names.txt"))

samples <- readr::read_table("sample_names.txt", col_names = F)%>%
  dplyr::pull(X1)

qfile_df <- data.table::fread(args[5])

# label Q file rownames and colnames
colnames(qfile_df) <- LETTERS[1:ncol(qfile_df)]

qfile_df$strain <- samples

qfile_long <- tidyr::gather(qfile_df, Pop, Perc, -strain) %>%
  dplyr::group_by(strain)%>%
  dplyr::filter(Perc == max(Perc))

POPSIZE <- length(unique(qfile_long$Pop))

POPGENOME_POPS <- list()
for(popgenome_pop in 1:POPSIZE){
  POPGENOME_POPS[[popgenome_pop]] <- as.character(dplyr::filter(qfile_long, Pop ==  LETTERS[popgenome_pop]) %>% dplyr::pull(strain))
}

#######################


GENOME_OBJECT <- PopGenome::readVCF(
  POPGENOME_VCF, 
  numcols = 1000, 
  tid = ANALYSIS_CHROM, 
  frompos = CHROM_START, 
  topos = CHROM_END, 
  gffpath = POPGENOME_GFF,
  include.unknown=TRUE)

GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)

GENOME_OBJECT <- PopGenome::set.populations(GENOME_OBJECT, POPGENOME_POPS, diploid = FALSE)

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
GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, keep.site.info = F, pi = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Fst Stats"))
GENOME_OBJECT <- PopGenome::F_ST.stats(GENOME_OBJECT, mode = "nucleotide", detail = TRUE)

save(GENOME_OBJECT, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_SUBPOPs_PopGenome_Object.Rda"))

# Process genome object

window_pos <- data.frame(get.neutrality(GENOME_OBJECT)[[1]]) %>%
  tibble::rownames_to_column() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = args[1],
                startWindow = as.numeric(strsplit(rowname,split = " ")[[1]][4]),
                endWindow = as.numeric(strsplit(rowname,split = " ")[[1]][6])) %>%
  dplyr::select(CHROM:endWindow)


pairFst <- data.frame(PopGenome::get.F_ST(GENOME_OBJECT, mode="nucleotide", pairwise = T)[[1]]) %>%
  dplyr::bind_cols(window_pos,.) %>%
  tidyr::gather(pops, Fst, -(CHROM:endWindow)) %>%
  tidyr::separate(pops, into = c("Pop1","Pop2"), sep = "\\.")

save(pairFst, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_SUBPOPs_FST_Statistics.Rda"))

system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

neutrality_ls <- list()
for(popnumber in 1:length(POPGENOME_POPS)){
  popname <- LETTERS[[popnumber]]
  
  n_df <- data.frame(get.neutrality(GENOME_OBJECT)[[popnumber]]) %>%
    dplyr::bind_cols(window_pos,.)
  
  d_df <- data.frame(get.diversity(GENOME_OBJECT)[[popnumber]]) %>%
    dplyr::bind_cols(window_pos,.)
  
  neutrality_ls[[popnumber]] <- dplyr::left_join(n_df,d_df, by = c("CHROM", "startWindow", "endWindow")) %>%
    dplyr::ungroup() %>%
    tidyr::gather(statistic, value, -(CHROM:endWindow)) %>%
    dplyr::mutate(Population = popname) %>%
    dplyr::select(CHROM:endWindow, Population, statistic, value)
  
}

neutrality_df <- dplyr::bind_rows(neutrality_ls)

save(neutrality_df, 
     file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_SUBPOPs_ND_Statistics.Rda"))
save(neutrality_df, 
     file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_SUBPOPs_ND_Statistics.Rda"))
