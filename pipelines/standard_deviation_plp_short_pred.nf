#!/usr/bin/env nextflow

/*
*  this module takes as input a txt file with the short ouptut prediction of phobius and
*  a corresponding plp file, if in the txt file are present more than one sequences the programm expects
*  the same axact sequence in the same order as the txt in the plp file, multiplp file
*  the multiplp file have to have this structure where every sequence prediction plp is separated by two line like this :
*  	# sp|Q8BI84|TGO1_MOUSE
*	#pos    aa      i       o       M       S       m       l
*  where the first comment line contains the id present in the line of the txt, remember the order has to be the same in the two files
*  this module computes the standard deviation value associated to what ever region of the prediciton is inputed
*  this region can be selected simply giving in the command line the  --KEYWORD argument
*  check the help message for the allowed values
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      this module takes as input a txt file with the short ouptut prediction of phobius option --INPUT_TXT and'
        log.info '      a corresponding plp file option --INPUT_PLP, if in the txt file are present more than one sequences the programm expects'
        log.info '      the same axact sequence in the same order as the txt in the plp file, multiplp file'
        log.info '      the multiplp file have to have this structure where every sequence prediction plp is separated by two line like this :'
        log.info '      	# sp|Q8BI84|TGO1_MOUSE'
        log.info '      	#pos    aa      i       o       M       S       m       l'
        log.info '      where the first comment line contains the id present in the line of the txt, remember the order has to be the same in the two files'
        log.info '      this module computes the standard deviation value associated to what ever region of the prediciton is inputed'
        log.info '      this region can be selected simply giving in the command line the  --KEYWORD argument'
        log.info '	the ouput format will be the following: '
        log.info '	average      standard deviation      id      range of segment    max'
        log.info '\n'
        exit 1
}



params.CONTAINER = "alessiovignoli3/tango-project:python_field_retr@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"
params.INPUT_TXT = "${params.TEST_DIR}bubbabubba"
params.INPUT_PLP = "${params.TEST_DIR}bubbabubba"
params.KEYWORD = false
params.MAX_ITER = false
params.SWITCH = "multi"					// flag used for the switch from search on a multi plp file or regular single plp
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.OUTPUT_NAME = "bubba"


// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"
include { ch_pairer } from "${params.PIPES}channel_files_pairer"
include { based_on_short_average_plp } from "${params.PIPES}average_plp_short_pred" 
include { average_unifier } from "${params.PIPES}average_plp_short_pred" addParams(SUFFIX: "stddeviation")



process based_on_short_stddeviation_plp {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if (params.MAX_ITER == false) filename
										else  null
										})
	tag { "${average_outfile}" }
	container params.CONTAINER

	input:
	tuple val(id), path(plp), path(average_outfile)
	path pyscript1
	val field_id

	output:
	path "${output_name}", emit: stddev_plp
	stdout emit: standardout

	script:
	prefix_plp = "${average_outfile}".split('\\.')[0]
	output_name = prefix_plp + ".stddeviation"
	//correct_plp = ""
	//for (item1 in plp) {
	//	if ( "${item1}".contains("${prefix_plp}") ) {
	//		correct_plp = item1
	//	}
	//}
	"""
	./${pyscript1} ${params.SWITCH} ${plp} ${field_id} ${average_outfile} ${params.MAX_ITER} 2>tmp.err 1>${output_name}
	cat tmp.err 
	"""
}


workflow short_pred_stddeviation_plp {
	
	take:
	input_txt
	input_plp
	field_keyword

	main:
	plp_average_computer_py = params.SCRIPTS + "plp_average_computer.py"
	stddeviation_computer_py = params.SCRIPTS + "plp_stddeviation_computer.py"
	pairer(input_txt, input_plp)
	based_on_short_average_plp(pairer.out.right_pairs, plp_average_computer_py, field_keyword)
	averageplp = based_on_short_average_plp.out.average_plp
	//println(in_txt.getClass())
	finalaverage = "there is nothing here, the --MAX_ITER has not been given a values, hence no need for this block to be executed this is not an error"
	in_plp = Channel.fromPath(input_plp)
	ch_pairer(in_plp, averageplp)
	//ch_pairer.out.view()
	based_on_short_stddeviation_plp(ch_pairer.out.correct_pairs, stddeviation_computer_py, field_keyword)
	stdeviationplp = based_on_short_stddeviation_plp.out.stddev_plp
	if ( params.MAX_ITER != false) {
		average_unifier(stdeviationplp.collect(), stddeviation_computer_py, field_keyword)
		finalaverage = average_unifier.out.final_average
	}

	emit:
	averageplp
	finalaverage
	stout1 = based_on_short_average_plp.out.standardout
	stout2 = based_on_short_stddeviation_plp.out.standardout
}



workflow {
	short_pred_stddeviation_plp(params.INPUT_TXT, params.INPUT_PLP, params.KEYWORD)
	//short_pred_stddeviation_plp.out.stout1.view()
	short_pred_stddeviation_plp.out.stout2.view()
	short_pred_stddeviation_plp.out.finalaverage.view()
	//short_pred_stddeviation_plp.out.averageplp.view()	
}
