# CePopulationGenetics-nf
archive of Stefan Z. CePopulationGenetics-nf repo with Tim C. edits

### Script versions
Two nextflow scripts are present in the repo.
- **cepopgen-nf.nf** is the original script.
- **cepopgen_TAC.nf** is the edited script. 
    - The pruned VCF used to perform the PCA is now output into `EIGESTRAT/INPUTFILES/PCA.vcf.gz`
    - The `eigenstrat_input.pedsnp` file now contains the correct snps.

### Usage notes
As of 20210801 only the PCA profile is recommended. To compare version outputs, run command with different `.nf` scripts.
```
nextflow run cepopgen_TAC.nf --out_base=163_PCA_0.1_hawaii_TAC --anc=XZ1516 --snv_vcf=../../Tim/Hawaii_popgen/20210120_Tim_Hawaii/popgen-20210120/WI.20210121.hard-filter.ref_strain_SNPs_only_fixup.vcf.gz -profile pca --pops=input_files/WI_20210121_163_Hawaii.tsv --eigen_ld=0.1
```
