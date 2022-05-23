#!/usr/bin/env nextflow

/*
*  this module takes as input a txt file with the short ouptut prediction of phobius and
*  a corresponding plp file, if in the txt file are present more than one sequences the programm expects
*  the same axact sequence in the same order as the txt in the plp file, multiplp file
*  the multiplp file have to have this structure where every sequence prediction plp is separated by two line like this :
*  	# sp|Q8BI84|TGO1_MOUSE
*	#pos    aa      i       o       M       S       m       l
*  where the first comment line contains the id present in the line of the txt, remember the order has to be the same in the two files
*  this module computes the average and the maximum plp value associated to what ever region of the prediciton is inputed
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
        log.info '      this module computes the average and the maximum plp value associated to what ever region of the prediciton is inputed'
        log.info '      this region can be selected simply giving in the command line the  --KEYWORD argument'
        log.info '	and it delimits the segment region to take into consideration it is in fct the main keyword label'
        log.info '	usefull when we want to know other averages for a given segment that are not the predicted label like,'
        log.info '	for example the selected main field keyword is s (special helix) but we want to compute also the average for'
        log.info '	the normal helix colum and loop for that predicted special helix segment'
        log.info '	use the --SEC_KEYWORD option, as for example  --KEYWORD s --SEC_KEYWORD n,c,l'
        log.info '	or --KEYWORD s --SEC_KEYWORD n,  this computes the average over the column in the plp for normal helix for the'
        log.info '	special helix segment that has beeen predicted on top of the average for such segment for special helix' 
        log.info '\n'
        exit 1
}



params.CONTAINER = "alessiovignoli3/tango-project@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b" // field retriever
params.INPUT_TXT = "${params.TEST_DIR}bubbabubba"
params.INPUT_PLP = "${params.TEST_DIR}bubbabubba"
params.KEYWORD = false
params.SEC_KEYWORD = false	// used when we want to know average on n for predicted s segment or any other possible combination of labels 
params.MAX_ITER = false
params.SWITCH = "multi"					// flag used for the switch from search on a multi plp file or regular single plp
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.SUFFIX = "average"


// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"



process based_on_short_average_plp {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if (params.MAX_ITER == false) filename
										else  null
										})
	tag { "${plp}" }
	container params.CONTAINER

	input:
	tuple val(id), path(txt), path(plp)
	path pyscript1
	val field_id

	output:
	path "${output_name}", emit: average_plp
	stdout emit: standardout

	script:
	prefix_plp = "${plp}".split('\\.')[0]
	output_name = prefix_plp + ".average"
	"""
	./${pyscript1} ${params.SWITCH} ${plp} ${field_id} ${txt} ${params.MAX_ITER} ${params.SEC_KEYWORD} 2>tmp.err 1>${output_name}
	cat tmp.err
	"""
}


process  average_unifier {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if ("${params.KEYWORD}" == "all") null
										else if ("${params.KEYWORD}".contains(',')) null
										else filename
										})
	tag { "${outname}" }
	container "ubuntu:latest@sha256:86ac87f73641c920fb42cc9612d4fb57b5626b56ea2a19b894d0673fd5b4f2e9"

	input:
	path list_of_averages
	path pyscript1
	val prefix
	
	output:
	path "${outname}", emit: final_average 
	stdout emit: standardout
	
	script:
	outname = "${params.INPUT_PLP}".split('\\.')[0].split('\\*')[0].split('/')[-1] + "-${prefix}.${params.MAX_ITER}tot_${params.SUFFIX}"
	"""
	for i in `echo ${list_of_averages}`; do cat `echo \$i | cut -d '[' -f 2 | cut -d ',' -f 1 | cut -d ']' -f 1` >> TMP ; done
	sort TMP | head -n ${params.MAX_ITER} > ${outname}
	#for i in `echo ${list_of_averages}`; do head `echo \$i | cut -d '[' -f 2 | cut -d ',' -f 1 | cut -d ']' -f 1`; done
	"""
}



workflow short_pred_average_plp {
	
	take:
	input_txt
	input_plp
	field_keyword

	main:
	//phobius_short_pred_field_retrieve_py = params.SCRIPTS + "phobius_short_prediction_field_retriever.py"
	plp_average_computer_py = params.SCRIPTS + "plp_average_computer.py"
	pairer(input_txt, input_plp)
	//pairer.out.right_pairs.view()
	based_on_short_average_plp(pairer.out.right_pairs, plp_average_computer_py, field_keyword)
	averageplp = based_on_short_average_plp.out.average_plp
	//println(in_txt.getClass())
	finalaverage = "there is nothing here, the --MAX_ITER has not been given a values, hence no need for this block to be executed this is not an error"
	if ( params.MAX_ITER != false) {
		average_unifier(averageplp.collect(), plp_average_computer_py, field_keyword)
		finalaverage = average_unifier.out.final_average
	}

	emit:
	averageplp
	finalaverage
	stout1 = based_on_short_average_plp.out.standardout
	//stout2 = average_unifier.out.standardout

}



workflow {
	short_pred_average_plp(params.INPUT_TXT, params.INPUT_PLP, params.KEYWORD)
	short_pred_average_plp.out.stout1.view()
	//short_pred_average_plp.out.stout2.view()
	short_pred_average_plp.out.finalaverage.view()
	short_pred_average_plp.out.averageplp.view()	
}
