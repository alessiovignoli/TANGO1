#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
		log.info ''
		log.info ''
		log.info ''
        log.info ''
		log.info ''
        log.info ''
		log.info ''
        log.info ''
		log.info ''
        log.info ''
        log.info '\n'
        exit 1
}



params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.OUTPUT_NAME = "bubba"
params.PUBLISH = true


process fasta_oneliner {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if (params.PUBLISH == true) filename
										else null
										})	
	tag { "${fasta}" }
	container params.CONTAINER

	input:
	tuple path(fasta), val(outname)
	path py_script1

	output:
	path "${outname}", emit: oneline_fasta
	//stdout emit: standardout

	script:
	"""
	python3 ${py_script1} ${fasta} ${outname}
	"""
}


workflow oneliner {
	
	take:
	pattern_to_fasta
	out_filename

	main:
	in_fasta = Channel.fromPath(pattern_to_fasta)
	list_in_and_out_names = in_fasta.map { [it, (it.baseName + ".${out_filename}")] }
	one_line_py = params.SCRIPTS + "one_line_per_fasta.py"
	fasta_oneliner(list_in_and_out_names, one_line_py)

	emit:
	onelinefasta = fasta_oneliner.out.oneline_fasta
	//stout = fasta_oneliner.out.standardout
}	


workflow {
	oneliner(params.INPUT, params.OUTPUT_NAME)
	//oneliner.out.onelinefasta.view()
	//oneliner.out.stout.view()
}


