#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
	log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
	log.info ""
        log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT_DOMAIN_INFO = "${params.TEST_DIR}bubbabubba.domain_info"
params.INPUT_TXT = "${params.TEST_DIR}bubbabubba.txt"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.FIELD = "s"

// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"

process  adjacency_finder {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${in_pred}" }

	input:
	tuple val(id), path(in_domain_info), path(in_pred)
	path py_script

	output:
	//stdout emit: standardout
	path "${out_name}_left.adj", emit: intermidiate_left_tmp
	path "${out_name}_center.adj", emit: intermidiate_center_tmp
	path "${out_name}_right.adj", emit: intermidiate_right_tmp

	script:
	out_name = "${in_domain_info}".split('\\.')[0]
	"""
	./${py_script} ${in_domain_info} ${in_pred} ${params.FIELD} ${out_name}
	"""
}


process summarizer_of_files {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
        container params.CONTAINER
        tag { "${all_one_type_annotations}" }

	input:
	path all_one_type_annotations
	
	output:
	stdout emit: standardout

	script:
	"""
	echo ${all_one_type_annotations}
	cat ${all_one_type_annotations} | cut -d ' ' -f5- | sort -u
	"""
}



workflow adjacent_domains_compiler {

	take:
	pattern_to_domaininfo
	patter_to_shortpred

	main:
	pairer(pattern_to_domaininfo, patter_to_shortpred)
	adj_finder_pyscript = params.SCRIPTS + "adjacent_annotation_finder.py"
	adjacency_finder(pairer.out.right_pairs, adj_finder_pyscript)
	adjacency_finder.out.intermidiate_left_tmp.collectFile(name: "sample_left.annot" ).set{ tmp_left }
	adjacency_finder.out.intermidiate_center_tmp.collectFile(name: "sample_center.annot" ).set{ tmp_center }
	adjacency_finder.out.intermidiate_right_tmp.collectFile(name: "sample_right.annot" ).set{ tmp_right }
	tmp_left.concat( tmp_center, tmp_right ).set{ proper_tmp }
	summarizer_of_files(proper_tmp)

	emit:
	stout = summarizer_of_files.out.standardout 
}

workflow {
	adjacent_domains_compiler(params.INPUT_DOMAIN_INFO, params.INPUT_TXT)
	adjacent_domains_compiler.out.stout.view()
}

