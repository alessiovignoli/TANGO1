#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info ''
        log.info ''
        log.info ''
        log.info '\n'
        exit 1
}



params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.PUBLISH = true


process fasta_volume_computer {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true, saveAs: { filename -> if (params.PUBLISH == true) filename
										else null
										})	
	tag { "${fasta}" }
	container params.CONTAINER

	input:
	tuple path(fasta), val(outname)
	path py_script1

	output:
	path "${outname}", emit: volume_file
	stdout emit: standardout

	script:
	"""
	python3 ${py_script1} -f ${fasta} -o ${outname}
	"""
}


workflow volume_computer {
	
	take:
	pattern_to_fasta

	main:
	in_fasta = Channel.fromPath(pattern_to_fasta)
	list_in_and_out_names = in_fasta.map { [it, (it.baseName + ".volumes")] }
	volume_computer_py = params.SCRIPTS + "fasta_volume_computer.py"
	fasta_volume_computer(list_in_and_out_names, volume_computer_py)

	emit:
	outfile = fasta_volume_computer.out.volume_file
	stout = fasta_volume_computer.out.standardout
}	


workflow {
	volume_computer(params.INPUT)
	volume_computer.out.outfile.view()
	volume_computer.out.stout.view()
}


