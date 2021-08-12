#!/usr/bin/env Rscript

# arguments
# 1 - Chromosome
# 2 - Ancestor: ${params.anc}
# 3 - VCF: ${vcf}
# 4 - GFF: ${params.gff}
# 5 - Strain population text file
# 6 - K size

# testing
# args <- c("I", "XZ1516", "/projects/b1059/projects/Stefan/CePopGen-nf/input_files/330_TEST.vcf.gz", "/projects/b1059/projects/Stefan/CePopGen-nf/input_files/Popgenome_files/gene_merge.gff", "/projects/b1059/projects/Stefan/CePopGen-nf/20181112-GATK4/ADMIXTURE/BEST_K/LD_0.1_MAF_0.004_BestK.Q", "6")

args <- commandArgs(trailingOnly = TRUE)

require(PopGenome)
require(data.table)
require(tidyverse)
require(WhopGenome)
require(bigmemory)

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

system(glue::glue("echo PopGenome - Reading VCF file")) 
GENOME_OBJECT <- PopGenome::readVCF(
  POPGENOME_VCF, 
  numcols = 10000, 
  tid = ANALYSIS_CHROM, 
  frompos = CHROM_START, 
  topos = CHROM_END, 
  approx = F, 
  gffpath = POPGENOME_GFF)

GENOME_OBJECT <- PopGenome::set.populations(GENOME_OBJECT, POPGENOME_POPS, diploid = FALSE)

GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)

GENOME_OBJECT <- PopGenome::splitting.data(GENOME_OBJECT, subsites= "exon")

system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Detail Stats"))
GENOME_OBJECT <- PopGenome::detail.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Neutrality Stats"))
GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT, detail = TRUE)
# system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
# GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Fst Stats"))
GENOME_OBJECT <- PopGenome::F_ST.stats(GENOME_OBJECT, mode = "nucleotide", detail = TRUE)

save(GENOME_OBJECT, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_WHOLE_POPULATION_PopGenome_Object.Rda"))

# Process genome object
windowStarts <- data.frame(Gene_Start = as.numeric(strsplit(GENOME_OBJECT@region.names, split = " - ")[[1]][1]),
                           Gene_End = as.numeric(strsplit(GENOME_OBJECT@region.names, split = " - ")[[1]][2]))

pairFst <- data.frame(PopGenome::get.F_ST(GENOME_OBJECT, mode="nucleotide", pairwise = T)[[1]],
                      windowStarts) %>%
  tidyr::gather(pops, Fst, -Gene_Start, -Gene_End) %>%
  tidyr::separate(pops, into = c("Pop1","Pop2"), sep = "\\.")

save(pairFst, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_WHOLE_POPULATION_FST_Statistics.Rda"))

system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

neutrality_ls <- list()
for(popnumber in 1:length(POPGENOME_POPS)){
  popname <- LETTERS[[popnumber]]
  neutrality_ls[[popnumber]] <- data.frame(PopGenome::get.neutrality(GENOME_OBJECT, theta = T, stats = T)[[popnumber]]) %>%
    dplyr::mutate(Population = popname,
                  Gene_Start = windowStarts$Gene_Start,
                  Gene_End = windowStarts$Gene_End) 
}

neutrality_df <- dplyr::bind_rows(neutrality_ls)
neutrality_df <- neutrality_df[,colSums(is.na(neutrality_df)) < nrow(neutrality_df)]
neutrality_df <- tidyr::gather(neutrality_df, statistic, value, -Population, -Gene_Start, -Gene_End)

save(neutrality_df, 
     file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_K{args[6]}_WHOLE_POPULATION_ND_Statistics.Rda"))
