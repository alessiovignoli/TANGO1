#!/usr/bin/env nextflow

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline computes "
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

params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster python3.9.5
params.REF = "${params.TEST_DIR}bubbabu"
params.IN = "${params.TEST_DIR}bubba"
params.FASTA = false
params.OUTPUT_DIR = "${params.TEST_DIR}"

// Modules dependencie section

include { retrieve_fastas } from "${params.PIPES}aln_generation_colouring"

/*
process   {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	//container "alessiovignoli3/tango-project@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"  field retriever
	container params.CONTAINER
	tag { "${}" }

	input:

	output:
	stdout emit: standardout

	script:
	"""
	"""
}



workflow  {

	take:

	main:

	emit:
	stout = summarizer_of_files.out.final_out
}
*/

workflow {
	if ( params.FASTA != false ) {
		retrieve_fastas(params.REF, params.FASTA)
		retrieve_fastas.out.retrieved_file.view()
	}
}

