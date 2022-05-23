#!/usr/bin/env nextflow

// ### PLEASE CHECK THE COMMENTS AND COMMENTED LINES IN THE SCRIPTS or SHELLS BLOCKS ###


// tm = tranmembrane
// om = Original model
// M3 = mark three third model created of the phobius model
// ns = negatie set
// ps = positive set TANGO1 proteins not used for the training of the model
// fs = frequency setthe TANGO1 proteins used for the assestment of the frequency of aa this prot are excluded from the positive set
// pred = prediction


/*
*  this module takes as input a fasta-like file with sequences spread on one or more than one line
*  the only important thing is that the headers are on one line and starting with > symbol
*  it basically outputs to the number of ines in total like wc -l and the number of headers
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      this module takes as input a fasta-like file with sequences spread on one or more than one line'
        log.info '      the only important thing is that the headers are on one line and starting with > symbol'
        log.info '      it basically outputs the number of ines in total like wc -l and the number of headers'
        log.info '\n'
        exit 1
}



params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT = "${params.TEST_DIR}bubbabubba"



process count_of_lines_and_headers {
	tag { "${fastas}" }
	container params.CONTAINER

	input:
	path fastas
	path py_script1
	
	output:
	stdout emit: standardout

	script:
	"""
	./${py_script1} ${fastas}
	"""
}
//count_of_lines_and_headers_out1.view()




workflow counter {
	
	take:
	pattern_to_fasta
	
	main:
	tmp = pattern_to_fasta.getClass()
	//println(tmp)
	in_fasta = pattern_to_fasta
	if ("${tmp}".contains("java.lang.String")) {
		in_fasta = Channel.fromPath(pattern_to_fasta)
	}
	counter_py = params.SCRIPTS + "huge_multifasta_counter.py"
	count_of_lines_and_headers(in_fasta, counter_py)

	emit:
	count = count_of_lines_and_headers.out.standardout
}


workflow {
	counter(params.INPUT)
	//counter.out.count.view()
}


