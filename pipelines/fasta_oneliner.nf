#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      this module takes as input a fasta-like file with sequens spread on more than one line and the output suffix to apend the input filename'
        log.info '      the only important thing is that the headers are on one line and starting with > symbol'
        log.info '      it basically outputs to the specified filename the same input fasta file but with each sequence on one line'
        log.info '\n'
	log.info '--INPUT		mandatory flag for one or more fasta files'
	log.info '--OUTPUT_NAME		optional flag for the output filename suffix to append to the name, files can not be fully renamed'
	log.info '			due to the fact that when in input there are many is impossibles for the user to specify all the relative oufilenames'
	log.info '			default value is -oneline.fa, so if the inoput file is example.test.fasta output will be example-oneline.fa'
	log.info '--OUTPUT_DIR		optional flag that specify the directory where the files are published=outputted, default value params.TEST_DIR'
	log.info '			the value of such variable can be found in the nextflow.config'
        exit 1
}



params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.OUTPUT_NAME = "-oneline.fa"
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
	list_in_and_out_names = in_fasta.map { [it, (it.simpleName + "${out_filename}")] }
	one_line_py = params.SCRIPTS + "one_line_per_fasta.py"
	fasta_oneliner(list_in_and_out_names, one_line_py)

	emit:
	onelinefasta = fasta_oneliner.out.oneline_fasta
	//stout = fasta_oneliner.out.standardout
}	

workflow oneliner_ch {

	take:
	channel_fasta
	out_filename

	main:
	list_in_and_out_names = channel_fasta.map { [it, (it.simpleName + "${out_filename}")] }
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


