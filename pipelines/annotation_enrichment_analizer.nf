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
        log.info " also --COLUMN_INFO can be used to specify wich column to consider with python notation where 0 is the first column (default)"
        log.info " like  --COLUMN_INFO 3 eill be the fourth column."
        log.info "	##### 		WRANING 	####"
        log.info " two other flags are necessary --TOTAL_IN and --TOTAL_BACK this two flags are used to compute the values for the table." 
        log.info " given the table that we want to obtain:  "
        log.info ""
        log.info "		man	woman"
	log.info "studying	a	b	a+b"
        log.info "not-studying	c	d	c+d"
	log.info "				a+b+c+d=n"
        log.info ""
	log.info ""
        log.info ""
        log.info "--TOTAL_IN is interpreted as a+b value that has to be associated to every line in the --INPUT_SET input file"
        log.info " while the value on each line in the column specified is interpeted as   a   "
        log.info ""
        log.info "simililarly --TOTAL_BACK is interpreted as c+d that is associated with --TOTAL_BACK background file"
        log.info "wile the value on each line in the column specified is interpeted as   c   this time "
        log.info ""
        log.info "this all pipeline takes as input a two (or more) files with the specifics said before and computes the table, for example"
        log.info "the line <24410 region disordered> of the input-file with --TOTAL_IN 55190   will be exact match searched using 'region disordered' on the background-file"
        log.info "let's say that this line is found:  <26116 region disordered> and --TOTAL_BACK 55190 is given"
        log.info "the reulting table will look like this:"
        log.info ""
	log.info "			region disordered	not region disordered"
        log.info "in-file		24410			30780"
        log.info "background-file	26116			29074"
	log.info ""
        log.info "the intermidiate output file line will look like <24410 30780 26116 29074 region disordered>. if b and c are inverted there is no difference for Fisher exact test"
	log.info ""
        log.info "now when two or more files are specified in --INPUT_SET or --BACKGROUND or both, all the cartesian product couples are generated and computed the fisher"
	log.info "example: --INPUT_SET -> in1.txt  in2.txt;  --BACKGROUND -> ba3.txt ba4.sc ba.xtx;  couples in1.txt,ba3.txt in1.txt,ba4.sc ... in2.txt,ba.xtx"
        log.info "no in1.txt,in2.txt couple nor ba3.txt,ba4.sc ecc.."
	log.info "in the above case --TOTAL_IN if it is one value is applied to all the input files, same goes for --TOTAL_BACK"
        log.info 'it is possible to specify for both the above something like --INPUT_SET "in*.txt"  --BACKGROUND "ba*" --TOTAL_IN 10,12 --TOTAL_BACK 50'
	log.info "in this case 10 refers to  in1.txt  and 12 refers to  in2.txt while 50 refers to all ba3.txt ba4.sc ba.xtx"
        log.info ""
	log.info "it need to be comma separated and it is up to the user making sure that the number given are as many as the inputs and in  the correct order"
        log.info 'try ls "in*.txt" before lounching the pipeline for the order in which it is presented to the sript'
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

include { fisher_exact } from "${params.PIPES}fisher_exact_test" addParams(CONTAINER: "alessiovignoli3/tango-project@sha256:259fd20c120c0e7691433a3ed54223f9e8d5b9ab7400e3eab4f4d3202a8c2594")	//  python slim-buster 3.9.5 and scipy 1.8

process table_values_computer  {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	//container "alessiovignoli3/tango-project@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"  field retriever
	container params.CONTAINER
	tag { "${in_set.baseName} + ${background.baseName}" }

	input:
	tuple path(in_set), val(in_set_total), path(background), val(background_total)
	val selected_column
	path pyscript

	output:
	path "*.ftable", emit: table_out
	stdout emit: standardout

	script:
	"""
	./${pyscript} ${in_set} ${in_set_total} ${background} ${background_total} ${selected_column}
	"""
}



workflow table_preparer {

	take:
	pattern_interest_set
	pattern_backgrounds
	total_in
	total_back
	column_to_use

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
		list_totalin = channel.of("${total_in}".split(',')).toList()
		matched_couples1 = channel.fromPath(pattern_interest_set).toSortedList().concat(list_totalin).toList().transpose(remainder: true)
		//matched_couples1.view()
		list_totalback = channel.of("${total_back}".split(',')).toList()
		matched_couples2 = channel.fromPath(pattern_backgrounds).toSortedList().concat(list_totalback).toList().transpose(remainder: true)
		//matched_couples2.view()
		coupled_filesets =  matched_couples1.combine(matched_couples2)
	}
	table_prepearer_pyscript = params.SCRIPTS + "fisher_table_preperer.py"
	table_values_computer(coupled_filesets, column_to_use, table_prepearer_pyscript)
	fisher_exact(table_values_computer.out.table_out, "two-sided")				//maybe need to put the mode in a vriable

	emit:
	stout = table_values_computer.out.standardout
	stout2 = fisher_exact.out.stout
	final_out = fisher_exact.out.final_fish
}

workflow {
	table_preparer(params.INPUT_SET, params.BACKGROUND, params.TOTAL_IN, params.TOTAL_BACK, params.COLUMN_INFO)
	table_preparer.out.stout.view()
	table_preparer.out.stout2.view()
	table_preparer.out.final_out.view()
}

