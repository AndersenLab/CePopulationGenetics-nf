#!/usr/bin/env nextflow 

date = new Date().format( 'yyyyMMdd' )
nextflow.preview.dsl=2


params.help                   = null
params.config                 = null
params.out_base               = null
params.anc                    = null
// params.ref                    = params.reference
params.cpu                    = "4"
params.snv_vcf                = params.snv_vcf
params.sv_vcf                 = params.sv_vcf
// params.admix_k                = params.admix_k
// params.admix_ld               = params.admix_ld
params.pops                   = params.pops
params.eigen_ld               = null
params.out                    = "${date}-${params.out_base}"

// params.best_Ks                = null
// params.haplotype              = params.haplotype
// params.popgene_window         = params.popgene_window
// params.popgene_window_subpop  = params.popgene_window_subpop
// params.popgene_gene           = params.popgene_gene
// params.popgene_gene_subpop    = params.popgene_gene_subpop

// params.admixture_100          = params.admixture_100
// params.admixture_seed         = params.admixture_seed
// params.phylo                  = params.phylo
params.species                = "c_elegans" //make sure to change species for trop or briggsae


if (params.help) {
    log.info '''
    ╔═╗   ╔═╗╦  ╔═╗╔═╗╔═╗╔╗╔╔═╗  ╔═╗╔═╗╔═╗╔═╗╔═╗╔╗╔  ╔╗╔╔═╗
    ║     ║╣ ║  ║╣ ║ ╦╠═╣║║║╚═╗  ╠═╝║ ║╠═╝║ ╦║╣ ║║║  ║║║╠╣ 
    ╚═╝o  ╚═╝╩═╝╚═╝╚═╝╩ ╩╝╚╝╚═╝  ╩  ╚═╝╩  ╚═╝╚═╝╝╚╝  ╝╚╝╚  
    '''
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow cepopgen-nf.nf --out_base Analysis --anc XZ1516"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--out_base             String                Name of folder to output results"
    log.info "--anc                  String                Name of ancestor strain to use"
    log.info "--snv_vcf              FILE                 Location to the small variant VCF to use for analysis"
    log.info "--pops                 FILE                 Location to sample population file"    
    log.info "--eigen_ld             String               Provide one or more LD values to test ('0.2' or '0.2,0.4,0.6'"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info "--species              String               c_elegans, c_briggsae, or c_tropicalis"
    // log.info "--ref                  FILE                 Location of the reference genome to use"
    log.info "--config               FILE                 Location of a custom configuration file"    
    // log.info "--sv_vcf               FILE                 Location to the structural variant VCF to use for analysis"
    // log.info "--admix_ld             STRING               Comma separated string of LDs to test for admixture analysis"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    // log.info "--email                STRING               email address for job notifications"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    // log.info "RAxML-ng               v0.5.1b"
    log.info "VCFtools               v0.1.16"
    // log.info "ADMIXTURE              v1.3"
    log.info "vcfanno                v0.2.8"
    log.info "EIGENSOFT              v6.1.4"
    // log.info "VARISCAN               v2.0"
    log.info "--------------------------------------------------------"    
    exit 1
} else {


log.info '''
╔═╗   ╔═╗╦  ╔═╗╔═╗╔═╗╔╗╔╔═╗  ╔═╗╔═╗╔═╗╔═╗╔═╗╔╗╔  ╔╗╔╔═╗
║     ║╣ ║  ║╣ ║ ╦╠═╣║║║╚═╗  ╠═╝║ ║╠═╝║ ╦║╣ ║║║  ║║║╠╣ 
╚═╝o  ╚═╝╩═╝╚═╝╚═╝╩ ╩╝╚╝╚═╝  ╩  ╚═╝╩  ╚═╝╚═╝╝╚╝  ╝╚╝╚  
'''
log.info ""
log.info "Species                                 = ${params.species}"
//log.info "Reference                               = ${params.ref}"
log.info "Ancestor                                = ${params.anc}"
log.info "Small Variant VCF                       = ${params.snv_vcf}"
log.info "Population File                         = ${params.pops}"
log.info "Max Population                          = ${params.admix_k}"
// log.info "LD to test                              = ${params.admix_ld}"
// log.info "cpu                                     = ${params.cpu}"
log.info "Output folder                           = ${params.out}"
// log.info "GFF3 File                               = ${params.gff}"
// log.info "Run Haplotype                           = ${params.haplotype}"
// log.info "Run Popgen window                       = ${params.popgene_window}"
// log.info "Run Popgen window subpops               = ${params.popgene_window_subpop}"
// log.info "Run Popgen gene subpops                 = ${params.popgene_gene_subpop}"
// log.info "Run Admixture cv100                     = ${params.admixture_100}"
// log.info "Run Admixture replicates                = ${params.admixture_seed}"
// log.info "Run Phylo                               = ${params.phylo}"
log.info ""
}


workflow {

    // start with VCF
    small_vcf = Channel.fromPath(params.snv_vcf)
    small_index = Channel.fromPath("${params.snv_vcf}" + ".tbi")

    // initialize population 
    File pop_file = new File(params.pops);
    pop_strains = Channel.from(pop_file.collect { it.tokenize( ' ' ) })
                 .map { POP, MAF, SM -> [POP, MAF, SM] }

    // extract ancestor
    small_vcf
      .combine(small_index) | extract_ancestor_bed


    // annotate small vcf
    small_vcf
      .combine(small_index)
      .combine(extract_ancestor_bed.out)
      .combine(pop_strains) | annotate_small_vcf 

    ld_range = Channel.of("${params.eigen_ld}")
                  .splitCsv()
                  .flatMap { it }

    // make vcf for eigenstrat - use LD provided
    annotate_small_vcf.out
      .combine(ld_range) | vcf_to_eigstrat_files

    vcf_to_eigstrat_files.out
      .combine(Channel.fromPath(params.eigen_par_no_removal)) | run_eigenstrat_no_outlier_removal

    vcf_to_eigstrat_files.out
      .combine(Channel.fromPath(params.eigen_par_outlier_removal)) | run_eigenstrat_with_outlier_removal


}


/*
==================================
~ > *                        * < ~
~ ~ > *                    * < ~ ~
~ ~ ~ > *  ANNOTATE VCF  * < ~ ~ ~
~ ~ > *                    * < ~ ~
~ > *                        * < ~
==================================
*/


/*
------------ Extract ancestor strain from the VCF and make bed file for annotations 
*/

process extract_ancestor_bed {

    publishDir "${params.out}/ANNOTATE_VCF", mode: 'copy'

    cpus 1

    input:
      tuple file(vcf), file(vcfindex)

    output:
      tuple file("ANC.bed.gz"), file("ANC.bed.gz.tbi")

      """
        bcftools query --samples ${params.anc} -f '%CHROM\\t%POS\\t%END\\t[%TGT]\\n' ${vcf} |\\
        awk -F"/" '\$1=\$1' OFS="\\t" |\\
        awk '{print \$1, \$2 = \$2 - 1, \$3, \$4}' OFS="\\t" |\\
        bgzip > ANC.bed.gz

        tabix ANC.bed.gz
        echo "ANCESTOR DONE"
      """
}

/*
------------ Annotate small variant VCF  
*/

process annotate_small_vcf {

    publishDir "${params.out}/ANNOTATE_VCF", mode: 'copy'

    conda '/projects/b1059/software/conda_envs/vcffixup'

    cpus 1

    input:
      tuple file(vcf), file(vcfindex), file("ANC.bed.gz"), file("ANC.bed.gz.tbi"), val(pop), val(maf), val(sm)

    output:
      tuple file("Ce330_annotated.vcf.gz"), file("Ce330_annotated.vcf.gz.tbi")


      """
        # get vcfanno files
        cp ${workflow.projectDir}/input_files/annotations/${params.species}/* .
        cat ${params.vcfanno_config} | sed 's/species/${params.species}/' > anno_config.toml

        bcftools view -s ${sm} ${vcf} -Oz -o ${pop}_pop.vcf.gz 

        vcfanno anno_config.toml ${pop}_pop.vcf.gz |\\
        awk '\$0 ~ "#" || \$0 !~ "Masked" {print}' |\\
        vcffixup - |\\
        bcftools filter -i N_MISSING=0 -Oz -o Ce330_annotated.vcf.gz

        tabix -p vcf Ce330_annotated.vcf.gz
      """
}



/*
------------ Prune VCF, Generate PLINK files
*/

process vcf_to_ped {

  tag {"PRUNING VCF FOR ADMIXTURE"}

  publishDir "${params.out}/ADMIXTURE/PLINK/", mode: 'copy'

  conda '/projects/b1059/software/conda_envs/vcffixup'


  input:
    tuple file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld)

  output:
    tuple file("ce_norm.vcf.gz"), file("ce_norm.vcf.gz.tbi"), val(nSM), val(maf), val(samples), val(ld), file("*.map"), file("*.ped"), file("plink.prune.in")

    """
    bcftools norm -m +snps ${vcf} -Oz -o ce_norm.vcf.gz
    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --maf ${maf} --set-missing-var-ids @:# --indep-pairwise 50 10 ${ld} --allow-extra-chr 
    
    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --maf ${maf} --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out LD_${ld}_MAF_${maf} --allow-extra-chr 
    """

}



/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  Run PCA and DAPC  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Prepare files for EIGENSTRAT
*/

process vcf_to_eigstrat_files {

  tag {"PREPARE EIGENSTRAT FILES"}

  conda '/projects/b1059/software/conda_envs/vcffixup'

  publishDir "${params.out}/EIGESTRAT/LD_${test_ld}/INPUTFILES", mode: 'copy'

  input:
    tuple file(vcf), file(vcfindex), val("test_ld")

  output:
    tuple file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi"), val(test_ld)


    """

    bcftools view --regions I,II,III,IV,V,X ${vcf} |\\
    bcftools norm -m +snps -Oz -o ce_norm.vcf.gz

    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${test_ld} --allow-extra-chr 

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out eigenstrat_input --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l ce_norm.vcf.gz |\\
    sort > sorted_samples.txt 

    bcftools view -v snps -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz -Oz -o PCA.vcf.gz
    
    tabix -p vcf PCA.vcf.gz

    bcftools view -v snps -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz |\\
    bcftools query -f '%CHROM\\t%CHROM:%POS\\t%cM\\t%POS\\t%REF\\t%ALT\\n' |\\
    sed 's/^III/3/g' |\\
    sed 's/^II/2/g' |\\
    sed 's/^IV/4/g' |\\
    sed 's/^I/1/g' |\\
    sed 's/^V/5/g' > eigenstrat_input.pedsnp      

    cut -f-6 -d' ' eigenstrat_input.ped |\\
    awk '{print 1, \$2, \$3, \$3, \$5, 1}'  > eigenstrat_input.pedind

    echo "rerun"
    """

}


/*
------------ Run EIGENSTRAT without removing outlier strains
*/

process run_eigenstrat_no_outlier_removal {

  publishDir "${params.out}/EIGESTRAT/LD_${test_ld}/NO_REMOVAL/", mode: 'copy'

  conda '/projects/b1059/software/conda_envs/vcffixup'

  input:
    tuple file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi"), val(test_ld), file(eigenparameters)

  output:
    tuple file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), \
    file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv")


    """

    smartpca -p ${eigenparameters} > logfile_no_removal.txt

    sed -n -e '/Tracy/,\$p' logfile_no_removal.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_no_removal.tsv
    """

}

/*
------------ Run EIGENSTRAT with removing outlier strains
*/

process run_eigenstrat_with_outlier_removal {

  conda '/projects/b1059/software/conda_envs/vcffixup'

  publishDir "${params.out}/EIGESTRAT/LD_${test_ld}/OUTLIER_REMOVAL/", mode: 'copy'

  input:
    tuple file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), \
    file ("markers.txt"), file ("sorted_samples.txt"), file ("PCA.vcf.gz"), file ("PCA.vcf.gz.tbi"), val(test_ld), file(eigenparameters)

  output:
    set file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), \
    file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), file("TracyWidom_statistics_outlier_removal.tsv")

   
    """
    smartpca -p ${eigenparameters} > logfile_outlier.txt

    sed -n -e '/Tracy/,\$p' logfile_outlier.txt |\
    sed -e '/kurt/,\$d' |\
    awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
    sed -e "s/[[:space:]]\\+/ /g" |\
    sed 's/^ //g' |\
    awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
    awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_outlier_removal.tsv
    """

}


