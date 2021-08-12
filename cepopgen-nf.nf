#!/usr/bin/env nextflow 

date = new Date().format( 'yyyyMMdd' )


params.help                   = null
params.config                 = null
params.out_base               = null
params.anc                    = null
params.ref                    = params.reference
params.cpu                    = "4"
params.snv_vcf                = params.snv_vcf
params.sv_vcf                 = params.sv_vcf
params.admix_k                = params.admix_k
params.admix_ld               = params.admix_ld
params.pops                   = params.pops
params.out                    = "${date}-${params.out_base}"

params.best_Ks                = null
params.haplotype              = params.haplotype
params.popgene_window         = params.popgene_window
params.popgene_window_subpop  = params.popgene_window_subpop
params.popgene_gene           = params.popgene_gene
params.popgene_gene_subpop    = params.popgene_gene_subpop

params.admixture_100          = params.admixture_100
params.admixture_seed         = params.admixture_seed
params.phylo                  = params.phylo

// use this to test more than one ld. 
// format of admix_ld = "ld\n0.1\n0.2\n0.3"
//ld_range = Channel.from(params.admix_ld)
//  .splitCsv(header:true)
//  .map{ row -> row.ld }
// lds = params.admix_ld.replaceAll(/\n/, ",").replaceAll(/ld,/, "")


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
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info "--ref                  FILE                 Location of the reference genome to use"
    log.info "--config               FILE                 Location of a custom configuration file"    
    log.info "--snv_vcf              FILE                 Location to the small variant VCF to use for analysis"
    log.info "--sv_vcf               FILE                 Location to the structural variant VCF to use for analysis"
    log.info "--pops                 FILE                 Location to sample population file"    
    log.info "--admix_k              INTEGER              Maximum admixture population size to test"
    log.info "--admix_ld             STRING               Comma separated string of LDs to test for admixture analysis"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    log.info "--email                STRING               email address for job notifications"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    log.info "RAxML-ng               v0.5.1b"
    log.info "VCFtools               v0.1.16"
    log.info "ADMIXTURE              v1.3"
    log.info "vcfanno                v0.2.8"
    log.info "EIGENSOFT              v6.1.4"
    log.info "VARISCAN               v2.0"
    log.info "--------------------------------------------------------"    
    exit 1
} else {


log.info '''
╔═╗   ╔═╗╦  ╔═╗╔═╗╔═╗╔╗╔╔═╗  ╔═╗╔═╗╔═╗╔═╗╔═╗╔╗╔  ╔╗╔╔═╗
║     ║╣ ║  ║╣ ║ ╦╠═╣║║║╚═╗  ╠═╝║ ║╠═╝║ ╦║╣ ║║║  ║║║╠╣ 
╚═╝o  ╚═╝╩═╝╚═╝╚═╝╩ ╩╝╚╝╚═╝  ╩  ╚═╝╩  ╚═╝╚═╝╝╚╝  ╝╚╝╚  
'''
log.info ""
log.info "Reference                               = ${params.ref}"
log.info "Ancestor                                = ${params.anc}"
log.info "Small Variant VCF                       = ${params.snv_vcf}"
log.info "Structural Variant VCF                  = ${params.sv_vcf}"
log.info "Population File                         = ${params.pops}"
log.info "Max Population                          = ${params.admix_k}"
log.info "LD to test                              = ${params.admix_ld}"
log.info "cpu                                     = ${params.cpu}"
log.info "Output folder                           = ${params.out}"
log.info "GFF3 File                               = ${params.gff}"
log.info "Run Haplotype                           = ${params.haplotype}"
log.info "Run Popgen window                       = ${params.popgene_window}"
log.info "Run Popgen window subpops               = ${params.popgene_window_subpop}"
log.info "Run Popgen gene subpops                 = ${params.popgene_gene_subpop}"
log.info "Run Admixture cv100                     = ${params.admixture_100}"
log.info "Run Admixture replicates                = ${params.admixture_seed}"
log.info "Run Phylo                               = ${params.phylo}"
log.info ""
}

K = Channel.from(2..params.admix_k)


ld_range = Channel.from(params.admix_ld)

/*
~ ~ ~ > * Define Contigs 
*/
CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]
Channel.from(CONTIG_LIST)
       .into{contigs_popgenome_window;
             contigs_popgenome_gene;
             contigs_popgenome_gene_subpops;
             contigs_popgenome_window_subpops}

// initialize population 
File pop_file = new File(params.pops);

pop_strains = Channel.from(pop_file.collect { it.tokenize( ' ' ) })
             .map { POP, MAF, SM -> [POP, MAF, SM] }

pop_strains
  .into{pop_ld;
        pop_subset}

pop_ld
  .combine(ld_range)
  .into{subset_vcf_input;
        plink_input; 
        print_plink_input }

// initialize VCF channel and spread to different channels for analysis
small_vcf = Channel.fromPath(params.snv_vcf)
small_index = Channel.fromPath("${params.snv_vcf}" + ".tbi")


// move this below VCF annotation when we decide to switch analysis to using masked VCF.

small_vcf
  .spread(small_index)
  .into { smallvcf_ancestor;
          smallvcf_annotations;
          smallvcf_skip_anno;
          }

// initialize eigenstrat parameter file chanel
eigenstrat_noremoval = Channel.fromPath(params.eigen_par_no_removal)
eigenstrat_removal = Channel.fromPath(params.eigen_par_outlier_removal)

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
      set file(vcf), file(vcfindex) from smallvcf_ancestor

    output:
      set file("ANC.bed.gz"), file("ANC.bed.gz.tbi") into anncestor_bed

    when:
      params.anc

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

    cpus 1

    input:
      set file(vcf), file(vcfindex) from smallvcf_annotations
      set file("ANC.bed.gz"), file("ANC.bed.gz.tbi") from anncestor_bed
      set val(pop), val(maf), val(sm) from pop_subset

    output:
      set file("Ce330_annotated.vcf.gz"), file("Ce330_annotated.vcf.gz.tbi") into annotated_vcf

    when:
      params.anc

      """
        bcftools view -s ${sm} ${vcf} -Oz -o ${pop}_pop.vcf.gz 

        vcfanno ${params.vcfanno_config} ${pop}_pop.vcf.gz |\\
        awk '\$0 ~ "#" || \$0 !~ "Masked" {print}' |\\
        vcffixup - |\\
        bcftools filter -i N_MISSING=0 -Oz -o Ce330_annotated.vcf.gz

        tabix -p vcf Ce330_annotated.vcf.gz
      """
}

if(!params.anc) {
	smallvcf_skip_anno
		.into{annotated_vcf;
		      print_vcf}

	params.popgen_anc = params.popgen_anc
} else {
	params.popgen_anc = params.anc
}

annotated_vcf
  .into { annotated_vcf_phylo;
          annotated_vcf_haplotype;
          annotated_vcf_popgen;
          annotated_vcf_eigenstrat;
          smallvcf_admixture;
          smallvcf_phylo;
          smallvcf_haplotype;
          smallvcf_classic_popgen_window;
          smallvcf_classic_popgen_gene;
          smallvcf_classic_popgen_gene_admix_pops;
          smallvcf_classic_popgen_window_admix_pops;
          smallvcf_haplotype_popgen;
          smallvcf_popstructure;
          smallvcf_eigenstrat}

/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > *  PHYLOGENY ANALYSIS  * < ~ ~ ~
~ ~ > *                          * < ~ ~ 
~ > *                              * < ~
========================================
*/

process vcf2phylip {

  publishDir "${params.out}/PHYLOGENY/PHYLIP", mode: "copy"

  input:
    set file(vcf), file(vcfindex) from annotated_vcf_phylo

  output:
    set file(vcf), file(vcfindex), file("*.phy") into phylip_out

  when:
    params.phylo

    """
    samples=`bcftools query -l ${vcf} | wc -l`
    python2 `which vcf2phylip.py` -i ${vcf} -m \$samples
    """
}

process generate_tree {

  cpus 8

  publishDir "${params.out}/PHYLOGENY/RAXML", mode: "copy"

  input:
    set file(vcf), file(vcfindex), file(phylip) from phylip_out

  output:
    set file("Population_Phylogeny.raxml.bestModel"),file("Population_Phylogeny.raxml.bestTree"),file("Population_Phylogeny.raxml.log") into raxml_tree

  when:
    params.phylo

    """
    raxml-ng --msa ${phylip} \\
      --model GTR \\
      --prefix Population_Phylogeny \\
      --threads ${task.cpus} 
    """
}

/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > *  ADMIXTURE ANALYSIS  * < ~ ~ ~
~ ~ > *                          * < ~ ~ 
~ > *                              * < ~
========================================
*/

smallvcf_admixture
  .spread(subset_vcf_input)
  .into { vcf_admix_plink;
          vcf_admix_print }


/*
------------ Prune VCF, Generate PLINK files
*/

process vcf_to_ped {

  tag {"PRUNING VCF FOR ADMIXTURE"}

  publishDir "${params.out}/ADMIXTURE/PLINK/", mode: 'copy'

  input:
    set file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld) from vcf_admix_plink

  output:
    set file("ce_norm.vcf.gz"), file("ce_norm.vcf.gz.tbi"), val(nSM), val(maf), val(samples), val(ld), file("*.map"), file("*.ped") into plink_output
    file("plink.prune.in") into pruned_marker_set

    """
    bcftools norm -m +snps ${vcf} -Oz -o ce_norm.vcf.gz
    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --maf ${maf} --set-missing-var-ids @:# --indep-pairwise 50 10 ${ld} --allow-extra-chr 
    
    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --maf ${maf} --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out LD_${ld}_MAF_${maf} --allow-extra-chr 
    """

}

/*
------------ Generate random seed for running ADMIXTURE 10 times
*/

Channel.from( 1000..2000 )
       .randomSample( 10, 234 )
       .into{ rseed_find_best_k;
              rseed_run_best_k}


/*
------------ Append population range and random seed to pruned plink files for ADMIXTURE analysis
*/

plink_output
  .into{ admixture_find_best_k;
         admixture_run_best_k}

K
  .spread(rseed_find_best_k)
  .spread(admixture_find_best_k)
  .into{ admixture_genome;
         admixture_regions}

/*
------------ Run ADMIXTURE analysis
*/

process run_admixture {

  tag { " ${pop} - ${rseed} "}

  publishDir "${params.out}/ADMIXTURE/${pop}/", mode: 'copy', pattern: '*.P'
  publishDir "${params.out}/ADMIXTURE/${pop}/", mode: 'copy', pattern: '*.Q'
  publishDir "${params.out}/ADMIXTURE/${pop}/", mode: 'copy', pattern: 'log_*'

  cpus 4

  input:
    set pop, rseed, file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped) from admixture_genome

  output:
    set pop, rseed, file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped), file("log_${pop}_${rseed}.out"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.P"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.Q") into admixture_output
    set pop, rseed, file("log_${pop}_${rseed}.out"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.P"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.Q") into admix_results

  when:
    params.admixture_seed

    """
      admixture --cv=10 -s ${rseed} ${ped} ${pop} -j4 | tee log_${pop}_${rseed}.out

      mv LD_${ld}_MAF_${maf}.${pop}.P LD_${ld}_MAF_${maf}_${pop}_${rseed}.P
      mv LD_${ld}_MAF_${maf}.${pop}.Q LD_${ld}_MAF_${maf}_${pop}_${rseed}.Q
    """
}

/*
------------ Combine 10 independent runs of each K
*/

admix_results
  .groupTuple()
  .into{ grouped_admix;
         grouped_print}

//grouped_print.println()

/*
------------ Combine log files of 10 independent runs of each K 
*/

process concat_replicate_logs {

  tag { " ${pop} "}

  echo true

  executor 'local'

  input:
    set pop, rseed, file(admixlog), file(admixp), file(admixq) from grouped_admix

  output:
    file("K${pop}_summary.txt") into concatenated_logs
  
  when:
    params.admixture_seed

  """
    grep -h CV log*.out |\\
    cut -f3- -d" " |\\
    sed 's/(\\|)\\|:\\|K=//g' > K${pop}_summary.txt
  """
}

/*
------------ Process CV results from ADMIXTURE analysis
*/

process concat_pop_logs {

  publishDir "${params.out}/ADMIXTURE/CV_Summary/", mode: 'copy'

  executor 'local'

  input:
    file(clog) from concatenated_logs.toSortedList()

  output:
    set file("admix_replicates_CV.tsv"), file("admix_summarized_CV.txt"), file("bestK.txt") into cv_summary
    file("bestK.txt") into bestk

  when:
    params.admixture_seed

  """
    # FULL RESULTS
    cat *summary.txt |\\
    sort -k1n |\\
    awk '\$1=\$1' OFS="\\t" |\\
    awk 'BEGIN{OFS="\\t"; print "K", "CV"}; {print \$0} OFS="\\t"' > admix_replicates_CV.tsv

    # Means of replicates
    cat *summary.txt |\\
    sort -k1n |\\
    awk '\$1=\$1' OFS="\\t" |\\
    awk 'BEGIN{OFS="\\t"; print "K", "CV"}; {print \$0} OFS="\\t"' |\\
    datamash -g 1 mean 2 -H |\\
    sed 's/GroupBy(\\|)\\|mean(//g' > admix_summarized_CV.txt

    # FIND BEST K - FIRST K WHERE NEXT HIGHER K CV VALUE IS SAME TO 2 decimal places
    cat *summary.txt |\\
    sort -k1n |\\
    awk '\$1=\$1' OFS="\\t" |\\
    awk 'BEGIN{OFS="\\t"; print "K", "CV"}; {print \$0} OFS="\\t"' |\\
    datamash -g 1 mean 2 -H |\\
    sed 's/GroupBy(\\|)\\|mean(//g' |\\
    awk 'NR>1{print \$0, sprintf("%3.2f", \$2-p)} {p = \$2}' |\\
    sed 's/-//g' |\\
    awk '\$3 == 0.00 {print}' |\\
    head -1 |\\
    cut -f-1 > bestK.txt

    # ADD PLUS MINUS 2 RANGE TO BEST K
    bk=`head -1 bestK.txt`

    START=\$(( 1+bk ))
    END=\$(( 2+bk ))
    for ((i=START;i<=END;i++)); do
        echo \$i >> bestK.txt
    done

    START=\$(( bk-2 ))
    END=\$(( bk-1 ))
    for ((i=START;i<=END;i++)); do
        echo \$i >> bestK.txt
    done
  """
}

/*
------------ Re-Run ADMIXTURE Analysis using Best K value 
*/

if(params.admixture_100){
  if (params.best_Ks) {

  File k_file = new File(params.best_Ks)
  best_k_file = Channel.from(k_file.collect { it.tokenize( ' ' ) })
                       .splitText(file: true)
                       

  admixture_run_best_k
    .spread(best_k_file)
    .into{ rerun_admixture;
           rerun_admixture_extra}

  } else {
    println """
    Please provide the population sizes you would like to run
    """
    System.exit(1)
  } 
} else {

    bestk
      .splitText(file: true)
      .set{k_range}

    admixture_run_best_k
      .spread(k_range)
      .into{ rerun_admixture;
             rerun_admixture_extra}
  }


process run_admixture_besk_k {

  publishDir "${params.out}/ADMIXTURE/BEST_K", mode: 'copy', pattern: '*.P'
  publishDir "${params.out}/ADMIXTURE/BEST_K", mode: 'copy', pattern: '*.Q'
  publishDir "${params.out}/ADMIXTURE/BEST_K", mode: 'copy', pattern: '*.ped'
  publishDir "${params.out}/ADMIXTURE/BEST_K", mode: 'copy', pattern: '*.map'

  cpus 4

  input:
    set file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped), file(k_size) from rerun_admixture

  output:
    set file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped), file(k_size), file("*.out"), file("*.P"), file("*.Q") into admixture_bestk
    file("*.Q") into admixture_bestk_to_popgenome

  when:
    params.admixture_100

    """
      bk=`head -1 ${k_size}`

      admixture --cv=100 ${ped} \$bk -j4 | tee log.\$bk.out
    """
}

admixture_bestk_to_popgenome
  .into{admixture_bestk_to_popgenome_gene;
        admixture_bestk_to_popgenome_window}

/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > *  HAPLOTYPE ANALYSIS  * < ~ ~ ~
~ ~ > *                          * < ~ ~
~ > *                              * < ~
========================================
*/

process_ibd=file("process_ibd.R")

/*
------------ Define IBDseq haplotype analysis parameters  
*/

minalleles = 0.01 // Specifies the minimum number of samples carrying the minor allele.
r2window = 1500 // Specifies the number of markers in the sliding window used to detect correlated markers.
ibdtrim = 0
r2max = 1

process ibdseq_haplotype {

  memory '64 GB' 

  publishDir "${params.out}/HAPLOTYPE", mode: 'copy'

  cpus 8

  input:
    set file(vcf), file(vindex) from smallvcf_haplotype

  output:
    file("haplotype.tsv") into haplotype_analysis

  when:
    params.haplotype

    """
      minalleles=\$(bcftools query --list-samples ${vcf} | wc -l | awk '{ print \$0*${minalleles} }' | awk '{printf("%d\\n", \$0+=\$0<0?0:0.9)}')
      if [[ \${minalleles} -lt 2 ]];
      then
          minalleles=2;
      fi;
      echo "minalleles=${minalleles}"
      for chrom in I II III IV V X; do
          java -jar `which ibdseq.r1206.jar` \\
              gt=${vcf} \\
              out=haplotype_\${chrom} \\
              ibdtrim=${ibdtrim} \\
              minalleles=\${minalleles} \\
              r2max=${r2max} \\
              nthreads=${task.cpus} \\
              chrom=\${chrom}
          done;

      cat *.ibd | awk '{ print \$0 "\\t${minalleles}\\t${ibdtrim}\\t${r2window}\\t${r2max}" }' > haplotype.tsv
    """
}

/*
------------ Stitch together IBDseq segments  
*/

process analyze_ibdseq {

    memory '64 GB' 

    publishDir "${params.out}/HAPLOTYPE", mode: 'copy'

    input:
        file("haplotype.tsv") from haplotype_analysis

    output:
        file("processed_haps.Rda")
        file("haplotype_plot_df.Rda") into plot_df

    when:
      params.haplotype

    """
      
      Rscript --vanilla `which process_ibd.R`
    """
}

/*
------------ Plot haplotypes, perform sweep analysis  
*/

process plot_ibdseq {

  memory '64 GB' 

  publishDir "${params.out}/HAPLOTYPE", mode: 'copy'

  input:
      file("haplotype_plot_df.Rda") from plot_df

  output:
      file("haplotype_length.pdf")
      file("max_haplotype_sorted_genome_wide.pdf")
      file("haplotype.pdf")
      file("sweep_summary.tsv")

  when:
    params.haplotype

  """
    Rscript --vanilla `which plot_ibd.R`
  """
}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  Classical Population Genetics Statistics  * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

process popgenome_whole_pop {

  memory '64 GB' 

  tag { CHROM }

  publishDir "${params.out}/POPGENOME/WINDOW/${CHROM}", mode: "copy"

  echo true

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_window
    each CHROM from contigs_popgenome_window

  output:
    set file("*Linkage_Statistics.Rda"), file("*Neutrality_Diversity_Statistics.Rda"), file("*WHOLE_POPULATION_Statistics.Rda") into popgenome_wholepop_statistics

  when:
    params.popgene_window

  script:
    """
      Rscript --vanilla `which PopGen_Window.R` ${CHROM} ${vcf} ${params.popgen_anc} ${params.popgenome_window} ${params.popgenome_slide} ${params.gff}
    """

}

process plot_popgenome_whole_pop {

  memory '64 GB' 

  publishDir "${params.out}/POPGENOME/PLOTS", mode: "copy", pattern: "*.png"
  publishDir "${params.out}/POPGENOME/WHOLE_GENOME", mode: "copy", pattern: "*.Rda"

  input:
    file("*") from popgenome_wholepop_statistics.collect()

  output:
    set file("Ce_Genome-wide_Neutrality_stats.Rda"), file("Ce_Genome-wide_Linkage_stats.Rda") into whole_genome_popgenome
    file("*.png") into whole_genome_popgenome_plots

  when:
    params.popgene_window

  script:
    """
      Rscript --vanilla `which Popgen_Plot_Window.R`
    """
}

process popgenome_window_kpop {

  tag { "${CHROM}" }

  memory '64 GB' 

  publishDir "${params.out}/POPGENOME/SUBPOPs/WINDOW/${CHROM}", mode: "copy"

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_window_admix_pops
    each file(popfile) from admixture_bestk_to_popgenome_window
    each CHROM from contigs_popgenome_window_subpops

  output:
    set val("${CHROM}"), file("*_FST_Statistics.Rda"), file("*_ND_Statistics.Rda") into fst_pops_window_output
    file("*_SUBPOPs_PopGenome_Object.Rda") into subpop_gobject

  when:
    params.popgene_window_subpop

    """
      ksize=`ls *.Q | rev | cut -d'.' -f 2 | rev`
      Rscript --vanilla `which PopGen_SubPops_Window.R` ${CHROM} ${params.popgen_anc} ${vcf} ${params.gff} ${popfile} \$ksize ${params.popgenome_window} ${params.popgenome_slide}
    """

}

process combine_popgenome_subpops {

  memory '64 GB' 

  publishDir "${params.out}/POPGENOME/SUBPOPs/WINDOW", mode: "copy", pattern: "*.Rda"

  input:
    file("*") from fst_pops_window_output.collect()

  output:
    set file("Ce_Genome-wide_Neutrality_stats.Rda"), file("Ce_Genome-wide_Fst_stats.Rda") into whole_genome_subpop_popgenome

  when:
    params.popgene_window_subpop

  script:

  """
    echo HELLO
    Rscript --vanilla `which PopGen_Concatenate_Subpops.R`
  """
}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  Classical PopGen Statistics - Gene Level * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

/*
------------ Calculate Gene-level PopGen Statistics
*/

process popgenome_gene_complete {

  memory '64 GB' 

  container 'andersenlab/popgen:v0.1'

  tag { "${CHROM}" }

  publishDir "${params.out}/POPGENOME/GENE/${CHROM}", mode: "copy"

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_gene
    each CHROM from contigs_popgenome_gene

  output:
    set val("${CHROM}"), file("CHROMOSOME-${CHROM}_WHOLE_POPULATION_Statistics.Rda"), file("CHROMOSOME-${CHROM}_WHOLE_POPULATION_PopGenome_Object.Rda") into fst_complete_gene_output

  when:
    params.popgene_gene

    """
      Rscript --vanilla `which PopGenome_Gene_exons.R` ${CHROM} ${params.popgen_anc} ${vcf} ${params.gff}
    """

}

process popgenome_gene_kpop {

  memory '64 GB' 

  tag { "${CHROM}" }

  publishDir "${params.out}/POPGENOME/SUBPOPs/GENE/${CHROM}", mode: "copy"

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_gene_admix_pops
    each file(popfile) from admixture_bestk_to_popgenome_gene
    each CHROM from contigs_popgenome_gene_subpops

  output:
    set val("${CHROM}"), file("*_WHOLE_POPULATION_PopGenome_Object.Rda") into fst_pops_gene_output

  when:
    params.popgene_gene_subpop

    """
      ksize=`ls *.Q | rev | cut -d'.' -f 2 | rev`
      Rscript --vanilla `which PopGenome_Gene_exons_subpops.R` ${CHROM} ${params.popgen_anc} ${vcf} ${params.gff} ${popfile} \$ksize
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

  publishDir "${params.out}/EIGESTRAT/INPUTFILES/", mode: 'copy'

  input:
    set file(vcf), file(vcfindex) from annotated_vcf_eigenstrat

  output:
    set file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind"), file("plink.prune.in"), file ("markers.txt") into eigenstrat_input

  when:
    params.eigenstrat_run

    """

    bcftools view --regions I,II,III,IV,V,X ${vcf} |\\
    bcftools norm -m +snps -Oz -o ce_norm.vcf.gz

    tabix -p vcf ce_norm.vcf.gz

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 ${params.eigen_ld} --allow-extra-chr 

    plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out eigenstrat_input --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l ce_norm.vcf.gz |\\
    sort > sorted_samples.txt 

    bcftools view -v snps -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz |\\
    bcftools query -f '%CHROM\\t%CHROM:%POS\\t%cM\\t%POS\\t%REF\\t%ALT\\n' ce_norm.vcf.gz |\\
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

eigenstrat_input
  .into{eigenstrat_no_outlier;
        eigenstrat_outlier_removal;
        eigenstrat_fst}


/*
------------ Run EIGENSTRAT without removing outlier strains
*/

process run_eigenstrat_no_outlier_removal {

  publishDir "${params.out}/EIGESTRAT/NO_REMOVAL/", mode: 'copy'

  input:
    set file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind") from eigenstrat_no_outlier
    file(eigenparameters) from eigenstrat_noremoval

  output:
    set file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv") into eigenstrat_outlier_removal_output

  when:
    params.eigenstrat_run

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


  publishDir "${params.out}/EIGESTRAT/OUTLIER_REMOVAL/", mode: 'copy'

  input:
    set file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind") from eigenstrat_outlier_removal
    file(eigenparameters) from eigenstrat_removal

  output:
    set file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), file("TracyWidom_statistics_outlier_removal.tsv") into eigenstrat_no_outlier_removal_output

  when:
    params.eigenstrat_run
   
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