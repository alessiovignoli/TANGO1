#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
	log.info ""
}


params.CONTAINER = 'lmansouri/phylo_3d:220225'
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.INPUT_ALN = false
params.INPUT_FASTA = false

params.gammaRate="1.0"
params.seedValue="5"
params.replicatesNum="100"
params.trimmer="trimal"
params.mode="8"


// include section 

include { align_generation; rename_gen } from "${params.PIPES}aln_generation_colouring" addParams(OUTPUT_DIR: "/dev/null")


process rename_aln {
	tag "${renamed_aln}"
	publishDir "${params.OUTPUT_DIR}", mode: 'copy', overwrite: false
	container "cbcrg/tcoffee@sha256:36c526915a898d5c15ede89bbc3854c0a66cef22db86285c244b89cad40fb855"

	input:
	path aln_file
        path rename_file
	
	output:
	path "${renamed_aln}", emit: renamed 
	stdout emit: standardout

	script:
	renamed_aln = "${rename_file}".split('\\.')[0] + '-renamed.aln'
	"""
	t_coffee -other_pg seq_reformat -in ${aln_file} -rename ${rename_file} -out ${renamed_aln}
	"""
}


process blocking_aln {
	tag"${id}"
	publishDir( "${params.OUTPUT_DIR}", mode: 'move', overwrite: false, saveAs: { filename -> if (filename.endsWith(".fa")) filename})
	container params.CONTAINER

	input:
	tuple val(id), path(in_aln)

	output:
	tuple val("${id}_trimmal_aln_auto"), path("*_auto.ph"), emit:  auto_trim
	tuple val("${id}_trimmal_aln_gappy"), path("*_gappy.ph"), emit: gappy_trim
	path "*.fa", emit: blocked_fasta
	stdout emit: standardout

	script:
	"""
	trimal -in ${in_aln} -out ${id}_trimmal_aln_auto.ph -phylip -automated1 --set_boundaries { 3247,3390 } -colnumbering 1>${id}_trimmal_auto.columnrelation
	trimal -in ${in_aln} -out ${id}_trimmal_aln_gappy.ph -phylip -gappyout -colnumbering 1>${id}_trimmal_gappy.columnrelation
	trimal -in ${in_aln} -out ${id}_trimmal_aln_auto.fa -fasta -automated1 
	trimal -in ${in_aln} -out ${id}_trimmal_aln_gappy.fa -fasta -gappyout
	trimal -in ${in_aln} -out ${id}_trimmal_aln.fa -fasta
	"""
}


process computing_ML_trees{
	errorStrategy 'ignore'
	tag"${id}"
	publishDir "${params.OUTPUT_DIR}", mode: 'copy', overwrite: false
	container params.CONTAINER

	input:
	tuple val(id), file(phylip) 

	output:
	tuple path("*.treefile"), path("*.boottrees"), emit: iqtree_ch
	stdout emit: standardout

	script:
	"""
	iqtree -s ${phylip} -b ${params.replicatesNum} 1>log
	"""
}


workflow auto_trim_ML {
	
	take:
	pattern_to_allignment
	patter_to_fasta

	main:
	input_aln = "" 
	if (params.INPUT_ALN != false) {
		input_aln = Channel.fromPath(pattern_to_allignment).map { item -> [ item.baseName.replace("_tmalign","") , item] }
	} else {
		input_fasta = Channel.fromPath(patter_to_fasta)
		align_generation(input_fasta)
		rename_gen(align_generation.out.aln_file)
		rename_aln(align_generation.out.aln_file, rename_gen.out.rename_file)
		input_aln = rename_aln.out.renamed.map { item -> [ item.baseName , item] }
	}
	blocking_aln(input_aln)	
	to_build_tree = blocking_aln.out.auto_trim.concat(blocking_aln.out.gappy_trim, input_aln) 
	computing_ML_trees(to_build_tree)

	emit:
	final_out = computing_ML_trees.out.iqtree_ch
	stout = computing_ML_trees.out.standardout
}


workflow {
	auto_trim_ML(params.INPUT_ALN, params.INPUT_FASTA)
	auto_trim_ML.out.final_out.view()
	auto_trim_ML.out.stout.view()
}
