#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline computes the fisher_exact_test using scipy library, it is mainly intended to be a module for other pipelines, but"
        log.info ""
        log.info "it is still launchable on it's own by terminal. The input files have to be like the following:"
        log.info "	6 2 1 4 proline rich domain"
        log.info "	43119 12071 38336 16854 transmembrane helical"
        log.info "a one or multiline file that has on one line the first 4 fields as integers values space separated and from the fifth on the name of the feature,"
        log.info "the feature being the annotation or the label associated with such values. The fopur values represent the table on which fisher"
	log.info "is done, first being topleft, second being top right, third being bottomleft and ecc.."
        log.info "it is responsability of the user to set them correctly, if it possible take a look at    annotation_enrichment_analizer.nf  "
        log.info "it prepares and lauches this  pipeline."
        log.info ""
        log.info "The input is given with --INPUT option, while the other optional flag is --NODE tat has as default value 'two-sided' "
        log.info "meaning  that is the probability a random table would have a probability equal to or less than the probability of the input table."
        log.info "this behavour can be changed giving as value 'greater' or 'less', where the first is the probability that a random table has x >= a,"
        log.info " which means with left top/ first number higher than what is given as input;"
        log.info "the seond instead is the one-sided p-value is the probability that a random table has x <= a."
        log.info "as many values can be specified as ones comma separated like: --MODE two-sided,greater,less"
        log.info "or any combination of the three."
        log.info ""
	log.info ""
        log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "alessiovignoli3/tango-project@sha256:259fd20c120c0e7691433a3ed54223f9e8d5b9ab7400e3eab4f4d3202a8c2594" // python slim-buster 3.9.5 and scipy 1.8
params.INPUT = "${params.TEST_DIR}bubbabu"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.MODE = "two-sided"

// Modules dependencie section

//include { pairer } from "${params.PIPES}input_files_pairer"


process fisher_exact_test  {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${input_table_file} + ${modes}" }
	//scratch true

	input:
	path input_table_file
	path pyscript
	each modes


	output:
	path "${out_name}", emit: fisher_out
	stdout emit: standardout

	script:
	out_name = ''
	if (modes == 'greater' || modes == 'less') {
		out_name = "${input_table_file.baseName}" + "_${modes}.fisher"
	} else {
		out_name = "${input_table_file.baseName}" + ".fisher"
	}
	"""
	./${pyscript} ${input_table_file} ${modes} 1>${out_name} 2>tmp
	cat tmp
	"""
}



workflow fisher_exact  {

	take:
	input_channel_files
	type_of_analysis

	main:
	list_of_modes = Channel.from("${type_of_analysis}".split(","))
	//list_of_modes.view()
	exact_fish_pyscript = params.SCRIPTS + "fisher_exact_test.py"
	fisher_exact_test(input_channel_files, exact_fish_pyscript,  list_of_modes)

	emit:
	final_fish = fisher_exact_test.out.fisher_out
	stout = fisher_exact_test.out.standardout
}

workflow {
	tmp_input = Channel.fromPath(params.INPUT)
	fisher_exact(tmp_input, params.MODE)
	fisher_exact.out.final_fish.view()
	fisher_exact.out.stout.view()
}

