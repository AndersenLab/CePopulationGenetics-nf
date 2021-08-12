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
GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, pi = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Linkage Stats"))
GENOME_OBJECT <- PopGenome::linkage.stats(GENOME_OBJECT, do.ZnS = TRUE, do.WALL = TRUE)

system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

save(GENOME_OBJECT, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Statistics.Rda"))

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Window Bins"))

windowStarts <- data.frame(snp_index = 1:length(colnames(GENOME_OBJECT@BIG.BIAL[[1]])),
                           position = colnames(GENOME_OBJECT@BIG.BIAL[[1]]))

slide_index <- cbind(data.frame(lapply(GENOME_OBJECT@SLIDE.POS,  function(x) as.numeric(floor(mean(x)))))) %>%
  tidyr::gather(temp, snp_index) %>%
  dplyr::select(-temp) %>%
  dplyr::left_join(., windowStarts, by = "snp_index")

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Linkage Stats"))

linkage_df <- data.frame(Wall.B = c(GENOME_OBJECT@Wall.B),
                         Wall.Q = c(GENOME_OBJECT@Wall.Q),
                         Kelly.Z_nS = c(GENOME_OBJECT@Kelly.Z_nS),
                         Rozas.ZZ = c(GENOME_OBJECT@Rozas.ZZ),
                         Rozas.ZA = c(GENOME_OBJECT@Rozas.ZA)) %>%
  dplyr::mutate(Population = "WHOLE_POPULATION",
                WindowPosition = slide_index$position) %>%
  tidyr::gather(LinkageStat, value, -Population, -WindowPosition)

save(linkage_df, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Linkage_Statistics.Rda"))

system(glue::glue("echo Generating PopGenome DataFrames - Extracting Neutrality and Diversity Stats"))

neutrality_df <- data.frame(Tajima.D = c(GENOME_OBJECT@Tajima.D),
                            n.segregating.sites = c(GENOME_OBJECT@n.segregating.sites),
                            Fu.Li.F = c(GENOME_OBJECT@Fu.Li.F),
                            Fu.Li.D = c(GENOME_OBJECT@Fu.Li.D),
                            Fay.Wu.H = c(GENOME_OBJECT@Fay.Wu.H),
                            Zeng.E = c(GENOME_OBJECT@Zeng.E),
                            theta_Tajima = c(GENOME_OBJECT@theta_Tajima),
                            theta_Watterson = c(GENOME_OBJECT@theta_Watterson),
                            theta_Achaz.Watterson = c(GENOME_OBJECT@theta_Achaz.Watterson),
                            theta_Achaz.Tajima = c(GENOME_OBJECT@theta_Achaz.Tajima),
                            theta_Fay.Wu = c(GENOME_OBJECT@theta_Fay.Wu),
                            theta_Zeng = c(GENOME_OBJECT@theta_Zeng),
                            nuc.diversity.within = c(GENOME_OBJECT@nuc.diversity.within),
                            PI = c(GENOME_OBJECT@Pi),
                            theta_Fu.Li = c(GENOME_OBJECT@theta_Fu.Li),
                            hap.diversity.within = c(GENOME_OBJECT@hap.diversity.within)) %>%
  dplyr::mutate(Population = "WHOLE_POPULATION",
                WindowPosition = slide_index$position) %>%
  tidyr::gather(statistic, value, -Population, -WindowPosition)

save(neutrality_df, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Neutrality_Diversity_Statistics.Rda"))