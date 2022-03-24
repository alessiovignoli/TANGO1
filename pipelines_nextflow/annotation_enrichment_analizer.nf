#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline computes a series of analizes the annotation frequency set comparing it to a background"
        log.info "one input, as many as you want backgrounds files, a necessary value that tell the total associated with each file, either one for all"
        log.info "or a list in the same order as the files, describe structure of files all should have same column aas info "
        log.info ""
        log.info "So the input set of annotation is passed with the flag --INPUT_SET whle the background with the flag --BACKGROUND"
        log.info "the structure of both files must be the following: a tab or space separated file with on each line the number and"
        log.info "the annotation name, for number is intended the real integer value of the interested annotation"
        log.info "something like:"
        log.info "	43119 transmembrane helical"
	log.info "	24410 region disordered"
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
        exit 1
}

params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster python3.9.5
params.INPUT_SET = "${params.TEST_DIR}bubbabu"
params.BACKGROUND =  "${params.TEST_DIR}bubbaback"
params.TOTAL_IN = 'bubba'
params.TOTAL_BACK = 'bub,'
params.COLUMN_INFO = 0		// python first value
params.OUTPUT_DIR = "${params.TEST_DIR}"

// Modules dependencie section

//include { pairer } from "${params.PIPES}input_files_pairer"


process table_values_computer  {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	//container "alessiovignoli3/tango-project@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"  field retriever
	container params.CONTAINER
	tag { "${}" }

	input:
	tuple path(in_set), val(in_set_total), path(background), val(background_total)

	output:
	stdout emit: standardout

	script:
	"""
	echo ${in_set} ${in_set_total} ${background} ${background_total}
	"""
}



workflow table_preparer {

	take:
	pattern_interest_set
	pattern_backgrounds
	total_in
	total_back

	main:
	coupled_filesets = ''
	list_totals = ("${total_in},${total_back}").split(',')
	//println(list_totals)
	if (list_totals.size() == 2) {
		int_set = channel.fromPath(pattern_interest_set).map{it -> [it, total_in]}
	        backgr = channel.fromPath(pattern_backgrounds).map{it -> [it, total_back]}
		coupled_filesets = int_set.combine(backgr)
	} else  {
		count1 = channel.fromPath(pattern_interest_set).count()
		count2 = channel.fromPath(pattern_backgrounds).count()
		num_files = count1.mix(count2).sum()
		//num_files.view()
		num_totals = channel.of(list_totals).count()
		//num_totals.view()
		tmp = num_files.mix(num_totals).collect()
		//tmp.view()
		tmp.map { if(it[0] != it[1]){ exit 1, "the number of input files and the number of values given are not the same, is up to the user to give them in the correct order and ammount, for more details pass --help as flag"
			 }}
		list_totalin = channel.of("${total_in}".split(',')).toSortedList()
		matched_couples1 = channel.fromPath(pattern_interest_set).toSortedList().mix(list_totalin).toSortedList().transpose(remainder: true)
		//matched_couples1.view()
		list_totalback = channel.of("${total_back}".split(',')).toSortedList()
		matched_couples2 = channel.fromPath(pattern_backgrounds).toSortedList().mix(list_totalback).toSortedList().transpose(remainder: true)
		//matched_couples2.view()
		coupled_filesets = matched_couples1.combine(matched_couples2)
	}
	//table_values_computer(coupled_filesets)

	emit:
	stout = coupled_filesets  // table_values_computer.out.standardout
}

workflow {
	table_preparer(params.INPUT_SET, params.BACKGROUND, params.TOTAL_IN, params.TOTAL_BACK)
	table_preparer.out.stout.view()
}

