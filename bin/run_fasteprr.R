#!/usr/bin/env Rscript

library(tidyverse)
library(FastEPRR)

args <- commandArgs(trailingOnly = TRUE)

dir.create("step1")
dir.create("step2")
dir.create("step3")

base_dir <- getwd()

vcf_file <- glue::glue("{base_dir}/{args[1]}.vcf.gz")

print(vcf_file)

FastEPRR_VCF_step1(vcfFilePath = "I.vcf.gz", 
                   winLength = "20", 
                   srcOutputFilePath = "{base_dir}/step1/")

FastEPRR_VCF_step2(srcFolderPath="{base_dir}/step1/", 
                   DXOutputFolderPath="{base_dir}/step2/")

FastEPRR_VCF_step3(srcFolderPath="{base_dir}/step1/", 
                   DXFolderPath="{base_dir}/step2/", 
                   finalOutputFolderPath=glue::glue("{base_dir}/{args[1]}_recombination.txt"))