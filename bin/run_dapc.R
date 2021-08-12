#!/usr/bin/env Rscript

# arguments
# 1 - SNPs
# 2 - Number of dicriminants to analyze
# 3 - Number of CPUs for analysis

# testing parameters
# args <- c("dapc_input.raw", 3, 1)

################################################################################
# initialize parameters

# load packages
library(adegenet)
library(tidyverse)
library(ape)

# colors for plotting
ancestry.colours <- c("gold2", "plum4","darkorange1", "lightskyblue2", 
                      "firebrick","burlywood3","gray51", 
                      "springgreen4", "lightpink2","deepskyblue4", 
                      "mediumpurple4", "orange", "maroon", "yellow3", "brown4", 
                      "yellow4","sienna4", "chocolate", "gray19")

# script arguments
args <- commandArgs(trailingOnly = TRUE)

# load snp set
snp_set <- adegenet::read.PLINK(args[1])

# set number of discriminant functions to retain
n_da <- as.numeric(args[2])

# set number of CPUs for analysis
ncpus <- as.numeric(args[3])

################################################################################

################################################################################
# run PCA analysis
ce_pca <- adegenet::glPca(snp_set, 
                          nf = 150,
                          n.cores = ncpus)

# define color for each strain based on first 3 PCs
myCol <- adegenet::colorplot(ce_pca$scores, 
                             ce_pca$scores, 
                             axes = c(1:3), 
                             transp=TRUE, cex=2)
abline(h=0,v=0, col="grey")

strain_colors <- data.frame(strain = snp_set$ind.names,
                            color = myCol)

save(ce_pca, file = "Adegent_PC_Object.Rda")
################################################################################

################################################################################
# find optimal PCs for DA

# determine groups
grp <- adegenet::find.clusters(snp_set, 
                               n.pca = 150, 
                               n.clust = n_da,
                               n.cores = ncpus)

# run DAPC
ce_dapc <- adegenet::dapc(snp_set, 
                          grp$grp, 
                          n.da = n_da, 
                          n.pca = 150,
                          glPca = ce_pca,
                          n.cores = ncpus)

# find optimal number of PCs to use for DA
find_optimal_pcs <- adegenet::optim.a.score(ce_dapc, n.cores = ncpus)

optimal_pcs <- find_optimal_pcs$best
################################################################################

################################################################################
# re-run DA analysis with optimal PC retention

# determine groups
grp <- adegenet::find.clusters(snp_set, 
                               n.pca = optimal_pcs, 
                               n.clust = n_da,
                               n.cores = ncpus)
# run DA
ce_dapc <- adegenet::dapc(snp_set, 
                          grp$grp, 
                          n.pca = optimal_pcs, 
                          n.da = n_da,
                          n.cores = ncpus)
# save DA object
save(ce_dapc, file = "Adegent_DAPC_Object.Rda")

# extract SNP loadings and save
snp_loading_df <- data.frame(Marker = as.character(snp_set$loc.names),
                             ce_dapc$var.contr) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CHROM = strsplit(as.character(Marker), split = ":")[[1]][1],
                POS = strsplit(strsplit(as.character(Marker), split = ":")[[1]][2], split = "_")[[1]][1]) %>%
  dplyr::select(-Marker) %>%
  dplyr::select(CHROM, POS, dplyr::contains("LD"))

write.table(snp_loading_df,
            file = "SNP_LOADINGS.tsv",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# extract Strain membership probabilities
memb_prob <- data.frame(strain = row.names(ce_dapc$posterior),
                        ce_dapc$posterior)

write.table(memb_prob,
            file = "STRAIN_MEMBERSHIP_PROBABILITIES.tsv",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

################################################################################

################################################################################
# plot DA results

pdf('DAPC_SCATTER.pdf',width = 10, height = 10)
scatter(ce_dapc, 
        col = ancestry.colours, 
        scree.da=FALSE, 
        bg="white", 
        pch=20,
        leg=TRUE, 
        txt.leg=paste("Cluster",1:n_da), 
        clab=0)
dev.off()

pdf('DAPC_LOADING.pdf',width = 20, height = 10)
loadingplot(ce_dapc$var.contr, thres=1e-3)
dev.off()

################################################################################

################################################################################
# Cross validation of DAPC analysis

# cv_snps <- tab(test_snps, NA.method="mean")
# 
# xval <- xvalDapc(cv_snps, 
#                  grp$grp,
#                  n.pca.max = 300, 
#                  training.set = 0.9,
#                  result = "groupMean", 
#                  center = FALSE, 
#                  scale = FALSE,
#                  n.pca = NULL, 
#                  n.rep = 30, 
#                  xval.plot = TRUE)
# 
# xval$`Number of PCs Achieving Lowest MSE`

