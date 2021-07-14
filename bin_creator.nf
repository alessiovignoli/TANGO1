#!/usr/bin/env nextflow

/*
*  this module takes as input a space or tab separated column file where the column requested contains numbers, by default is the first
*  all other columns will not be considered
*  it basically outputs to stdout or a file if specified a python3 dyctionary where the keys will be the intervals and the argument the 
*  amount of numbers (rows) fell into that bin
*  the number of bins can be manually specified using the option --BIN_NUM, is set to ten by default
*  this module does require the range of the numbers, basicallly the minimun and the maximum by deafault isassumed to be between 0 and 1
*  this can be modified with the argument --RANGE num1,num2
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      this module takes as input a space or tab separated column file option --INPUT where the column requested contains numbers option --COLUMN'
        log.info '      by default is the first with value 0, it follows the python notation where the first is actually number 0 the second number 1 ecc.'
        log.info '	all other columns will not be considered'
        log.info '      it basically outputs to stdout or a file option --SWITCH file, if specified a python3 dyctionary where the keys will be the intervals and the argument the'
        log.info '	amount of numbers (rows) fell into that bin'
        log.info '	the number of bins can be manually specified using the option --BIN_NUM, is set to ten by default'
        log.info '	this module does require the range of the numbers, basicallly the minimun and the maximum by deafault isassumed to be between 0 and 1'
        log.info '	this can be modified with the argument --RANGE num1,num2'
        log.info '	num1 and num2 are expected to be either integers or python floats (4.5 for example no comma just dot)'
        log.info '	when launching from command line if a num has to be negative put \\ in front of the minus sign like \\-3.2,3.2'
        log.info '\n'
        exit 1
}



params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.COLUMN = 0
params.BIN_NUM = 10
params.RANGE = "0,1"
params.SWITCH = "stdout"
params.OUTPUT_DIR = "${params.TEST_DIR}"


process bin_compiler {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if (params.SWITCH == "stdout") null
										else  filename
										})
	tag { "${infile}" }
	container params.CONTAINER

	input:
	path infile
	path py_script1
	val column

	output:
	path "${output_name}", emit: dict_of_bins
	stdout emit: standardout

	script:
	output_name = "${infile}".split('\\.')[0] + "." + "${infile}".split('\\.')[1] + "_bins"
	if (params.SWITCH == "stdout") {
		"""
		touch ${output_name}
		./${py_script1} ${infile} ${column} ${params.BIN_NUM} ${params.RANGE} 2>tmp.err 
		cat tmp.err
		"""
	}
	else {
		"""
		./${py_script1} ${infile} ${column} ${params.BIN_NUM} ${params.RANGE} 2>tmp.err  1>${output_name}
		cat tmp.err
		"""
	}
}




workflow  biner {
	
	take:
	pattern_to_input
	desired_column

	main:
	in_file = Channel.fromPath(pattern_to_input)
	bin_computer_py = params.SCRIPTS + "bin_computer.py"
	bin_compiler(in_file, bin_computer_py, desired_column)
	//bin_uniter(bin_compiler.out.dict_of_bins.collect())

	emit:
	stout = bin_compiler.out.standardout
	//stout = bin_uniter.out.standardout
}	


workflow {
	biner(params.INPUT, params.COLUMN)
	//oneliner.out.onelinefasta.view()
	biner.out.stout.view()
}


