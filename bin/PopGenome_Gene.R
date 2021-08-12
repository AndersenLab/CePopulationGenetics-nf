#!/usr/bin/env Rscript

# arguments
# 1 - Chromosome
# 2 - Ancestor: ${params.anc}
# 3 - VCF: ${vcf}
# 4 - GFF: ${params.gff}
# 5 - 

# testing
# args <- c("I", "XZ1516", "/projects/b1059/projects/Stefan/CePopGen-nf/input_files/330_TEST.vcf.gz", "/projects/b1059/projects/Stefan/CePopGen-nf/input_files/Popgenome_files/gene_merge.gff")

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

system(glue::glue("echo PopGenome - Reading VCF file")) 

vcf_handle <- vcf_open(POPGENOME_VCF)

GENOME_OBJECT <- Whop_readVCF(
  vcf_handle, 
  numcols = 1000, 
  tid = ANALYSIS_CHROM, 
  frompos = CHROM_START, 
  topos = CHROM_END, 
  gffpath = POPGENOME_GFF)

GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)

GENOME_OBJECT <- PopGenome::splitting.data(GENOME_OBJECT, subsites= "gene")

system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Detail Stats"))
GENOME_OBJECT <- PopGenome::detail.stats(GENOME_OBJECT)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Neutrality Stats"))
GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT, detail = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, pi = TRUE)
system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Fst Stats"))
GENOME_OBJECT <- PopGenome::F_ST.stats(GENOME_OBJECT, mode = "nucleotide", detail = TRUE)
#system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Linkage Stats"))
#GENOME_OBJECT <- PopGenome::linkage.stats(GENOME_OBJECT, do.ZnS = TRUE, do.WALL = TRUE)

system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

save(GENOME_OBJECT, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_PopGenome_Object.Rda"))

# Object to data frame
gene_names <- list()

for(g in 1:length(GENOME_OBJECT@region.names)){
  gene_names[[g]] <- data.frame(Population = "Whole Population",
                     Outgroup = OUTGROUP,
                     CHROM = ANALYSIS_CHROM,
                     Gene_Start = as.numeric(strsplit(GENOME_OBJECT@region.names[g], split = " - ")[[1]][1]),
                     Gene_End = as.numeric(strsplit(GENOME_OBJECT@region.names[g], split = " - ")[[1]][2]),
                     n_sites = GENOME_OBJECT@n.segregating.sites[g],
                     Watterson_theta = GENOME_OBJECT@theta_Watterson[g],
                     Tajima_theta = GENOME_OBJECT@theta_Tajima[g],
                     Achaz.Watterson_theta = GENOME_OBJECT@theta_Achaz.Watterson[g],
                     Achaz.Tajima_theta = GENOME_OBJECT@theta_Achaz.Tajima[g],
                     Fay.Wu_theta = GENOME_OBJECT@theta_Fay.Wu[g],
                     Zeng_theta = GENOME_OBJECT@theta_Zeng[g],
                     Fu.Li_theta = GENOME_OBJECT@theta_Fu.Li[g],
                     TajimaD = GENOME_OBJECT@Tajima.D[g],
                     Fu.Li.F = GENOME_OBJECT@Fu.Li.F[g],
                     Fu.Li.D = GENOME_OBJECT@Fu.Li.D[g],
                     Fay.Wu.H = GENOME_OBJECT@Fay.Wu.H[g],
                     Zeng.E = GENOME_OBJECT@Zeng.E[g],
                     #Wall.B = GENOME_OBJECT@Wall.B[g],
                     #Wall.Q = GENOME_OBJECT@Wall.Q[g],
                     #Kelly.Z_nS = GENOME_OBJECT@Kelly.Z_nS[g],
                     #Rozas.ZZ = GENOME_OBJECT@Rozas.ZZ[g],
                     #Rozas.ZA = GENOME_OBJECT@Rozas.ZA[g],
                     hap.diversity.within = GENOME_OBJECT@hap.diversity.within[g],
                     Fst = as.numeric(GENOME_OBJECT@nucleotide.F_ST)[g]
  )
}

gene_stats_df <- dplyr::bind_rows(gene_names) %>%
  tidyr::gather(Statistic, Value, -(Population:Gene_End))

save(gene_stats_df, 
     file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Statistics.Rda"))
