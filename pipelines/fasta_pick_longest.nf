#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline picks the longest fasta sequence and corresponding header and saves it to a given file"
        log.info "in case of multiple input fasta files it picks the longest from each one of them and saves to one output file"
        log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info '--IN		input fasta or fasta filepath with glob pattern, like "../dir/some_file*.fasta"'
        log.info "		in case of many files the longest sequence is picked from each one"
        log.info "--OUT		mandatory flag, the name or path of the file (non existant) in which the selected sequences have to be written"
	log.info "		no need to add a .fasta or similar .fasta is already added"
        log.info "		if the path is a glob pattern pieline output is uncertain."
        log.info "		the path of the output file can be absolute, or relative to where the pipeline is lounched"
	log.info "		through the use of the variable of nextflow   launchDir   "
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster python3.9.5
params.IN = "${params.TEST_DIR}bubbabu"
params.OUT = false


// Modules dependencie section

include { converger } from "${params.PIPES}output_files_uniter" addParams(OUTPUT_DIR: "${launchDir}", SWITCH_PUBLISH: "bubba")


process fasta_longest_picker  {
	container params.CONTAINER
	tag { "${fasta}" }

	input:
	path fasta
	path pyscript

	output:
	path "*.tmp", emit: tmp_files
	//stdout emit: standardout	for debug

	script:
	"""
	python3 ${pyscript} --fasta ${fasta} > ${fasta}.tmp
	"""
}


workflow  fasta_longest {

	take:
	pattern_to_fastas
	pattern_to_out

	main:

	// error section of missing or wrong inputs

	if ( !params.OUT ) {
		log.info "no --OUT flag specified please give a vild path, or type --help for description of pipeline"
	}

	fastas = channel.fromPath(pattern_to_fastas)
	fasta_pick_the_longest_pyscript = params.SCRIPTS + "fasta_pick_the_longest.py"
	fasta_longest_picker(fastas, fasta_pick_the_longest_pyscript)
	converger(fasta_longest_picker.out.tmp_files, pattern_to_out, "fasta")

	emit:
	stout = converger.out.stout
}

workflow {
	fasta_longest(params.IN, params.OUT)
	fasta_longest.out.stout.view()
}

