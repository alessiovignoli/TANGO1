#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline retrieves the ids present in a phobius short output prediction based on a series ofreasoning"
        log.info "in brief it checks if in a prediction of a protein there is n number of consecutive feature and extract the id"
        log.info "the input is only a short phobius prediction redirected to a file --INPUT_TXT flag, it has the fllowing structure:"
        log.info "	UniRef50_A0A3D5HTR8             7  Y n3-13c17/18o41-n-60i217-n-237o257-n-277i298-n-324o330-n-351i363-n-383o389-n-409i"
        log.info "	UniRef50_UPI000CE2750B          3  0 i42-n-59o79-n-96i223-n-244o"
        log.info "the prediction section and features ranges are extracted using     phobius_short_prediction_field_retriever.py"
        log.info "then the distance between the requested feature is checked to be less than a threshold more later,"
        log.info "if each segment is below the theshold at distance to its following it is considered consecutive and if all n segments"
        log.info "respect this rule the id is selected and wrote to the output file"
	log.info "--FIELD is used to specify wich is the feature of the prediction to look for"
	log.info 'default is s for special helix, all the possible are: [c, i, o, -, n, s, l]  c = signal peptide, i = inside membrane(cytoplasm),'
        log.info 'o = outside membrane, - = helix (in phobius originalmodel), (only in phobius-M7or later) => -n- = normal-helix'
        log.info '-s- = special-helix and -l- = loop-inramembrane'
        log.info 'they have to be given to the pipeline in one letter code'
        log.info "the - has to be given like \\-"
        log.info ""
        log.info "--NUMBER_CONSEC is used as option to specify how many feature (n) are to look for, if 2 only protein with 2 normal conseciìutive helices are selected"
        log.info ""
        log.info "--THRESHOLD this is the length to be less or equal to be considered consecutive, default 25"
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

