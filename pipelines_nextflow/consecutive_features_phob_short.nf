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
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        exit 1
}

params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster
params.INPUT_TXT = "${params.TEST_DIR}bubbabubba.txt"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.FIELD = "s"
params.NUMBER_CONSEC = false	// default here would be 2
params.THRESHOLD = "25"		// default in the python script



process  consecutive_finder {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container "alessiovignoli3/tango-project@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b" // field retriever
	tag { "${in_pred}" }

	input:
	path in_pred
	val field_pred
	path py_script

	output:
	stdout emit: standardout

	script:
	out_name = "${in_pred}".split('\\.')[0] + ".ids"
	if (params.NUMBER_CONSEC != false) {
		"""
		./${py_script} ${in_pred} ${out_name} ${field_pred} ${params.NUMBER_CONSEC} ${params.THRESHOLD}
		"""
	} else {
		"""
		./${py_script} ${in_pred} ${out_name} ${field_pred}
		"""
	}
}


workflow exact_consecutive_id_retriever {

	take:
	patter_to_shortpred
	field_value

	main:
	in_shortpred = Channel.fromPath(patter_to_shortpred)
	consec_finder_pyscript = params.SCRIPTS + "consecutive_segments_id_retriever.py"
	consecutive_finder(in_shortpred, field_value, consec_finder_pyscript)

	emit:
	stout = consecutive_finder.out.standardout 
}

workflow {
	exact_consecutive_id_retriever(params.INPUT_TXT, params.FIELD)
	exact_consecutive_id_retriever.out.stout.view()
}

