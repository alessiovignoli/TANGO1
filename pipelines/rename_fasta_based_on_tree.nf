#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline order a multi fasta file on the basis of a tree in newick format "
        log.info "the fasta headers have to contain at least one univocal id that is used in the tree to indicate that leave"
        log.info "for example the newick trees indicates the leaves with the taxid number like human is represented as 9606"
        log.info "now the multifasta file has to have only one header with such id somwhere in it for it to be correctly ordered"
        log.info "same goes for all othere sequences, they need to have a one to ne relationship with the id found in the tree file"
        log.info "it is also possible to specify a --DELIMITER flag that specifies what is delimiting the field like:"
        log.info '--DELIMITER " "   or   --DELIMITER "OX="   quotations are mandatory you can use both single and double'
        log.info "the delimiter will be tried to use to isolate the univocal id found in the tree in the header of the fasta, two"
        log.info "combinations are tried <delimiterIDdelimiter> <delimiterID> like the first example < ID > < ID> or the second"
	log.info "<OX=IDOX=> <OX=ID>   for more detail read  fasta_renamer_based_on_tree.py"
        log.info ""
        log.info "with the flag --INPUT_TREE the tree is given to the pipeline, with --INPUT_FASTA the corresponding fasta file is given"
        log.info "for multiple couples of tree and fasta there is the following sheme:"
        log.info ""
	log.info "####  WARNIONG #####"
	log.info "This module imports channel_files_pairer.nf ch_pairer workflow to handle matching of multiple file pairs."
	log.info "go look help section of that script for more detaqils."
        log.info ""
        log.info "to work on the fasta file the sequence has to be on one line only so by default the pipeline will do it by itself"
        log.info "this behaviour can be changed, so to skip this process, pass any value to --ONE_LINE flag"
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster python3.9.5
params.INPUT_TREE = "${params.TEST_DIR}bubbabu"
params.INPUT_FASTA = "${params.TEST_DIR}bubb"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.ONE_LINE = false
params.DELIMITER = ""

// Modules dependencie section

//include { pairer } from "${params.PIPES}input_files_pairer"
include { oneliner } from "${params.PIPES}fasta_oneliner" addParams(PUBLISH: "false")
include { ch_pairer } from "${params.PIPES}channel_files_pairer"



process   matcher_renamer_fasta {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true)
	container params.CONTAINER
	tag { "${matcher}" }

	input:
	tuple val(matcher), path(treeFile), path(fastaFile)
	val delimit
	path matchrename_script

	output:
	path "${out_renamedfasta}", emit: renamed_fasta
	path "correspondance_fasta_tree.tab", emit: correspond_table
	stdout emit: standardout

	script:
	out_renamedfasta = "${fastaFile}".split('\\.')[0] + ".renamed.fasta"
	if (delimit.isEmpty()) {
		"""
		python3 ${matchrename_script} -t ${treeFile} -f ${fastaFile} -o ${out_renamedfasta}
		"""
	} else {
		"""
		python3 ${matchrename_script} -t ${treeFile} -f ${fastaFile} -o ${out_renamedfasta} -d ${delimit}
		"""
	}
}



workflow  match_and_order_treefasta {

	take:
	input_tree
	input_fasta
	delimiter

	main:
	// scripts
	matchrename_pyscript = params.SCRIPTS + "fasta_renamer_based_on_tree.py"
	

	in_fasta = ''
	if (params.ONE_LINE != false) {
                in_fasta = Channel.fromPath(input_fasta)
        } else {
                suffix = "oneline"
                in_fasta = oneliner(input_fasta, suffix)
	}
	in_tree = Channel.fromPath(input_tree)
	ch_pairer(in_tree, in_fasta)
	matcher_renamer_fasta(ch_pairer.out.correct_pairs, delimiter, matchrename_pyscript)

	emit:
	final_out = matcher_renamer_fasta.out.renamed_fasta
	stout = matcher_renamer_fasta.out.standardout
}

workflow {
	match_and_order_treefasta(params.INPUT_TREE, params.INPUT_FASTA, params.DELIMITER)
	match_and_order_treefasta.out.final_out.view()
	match_and_order_treefasta.out.stout.view()
	
}

