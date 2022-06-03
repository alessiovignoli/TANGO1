#!/usr/bin/env nextflow

/*
*  this module takes as input a fasta-like file and the number of entries to be selected 
*  if the input is instead a glob pattern for many fasta-like files the number of entries to be selected has to be the total one
*  the ammount of sequences to be extracted per each file is computed dinamically 
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
	log.info 'This is the help section of this pipeline'
	log.info 'Here is a brief description of the pipeline itself  : '
	log.info '	this module takes as input a fasta-like file --INPUT and the number of entries to be selected --NUMBER'
	log.info '	if the input is instead a glob pattern for many fasta-like files the number of entries to be selected has to be the total one'
	log.info '	the ammount of sequences to be extracted per each file is computed dinamically'
	log.info '\n'
	exit 1
}



params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT = "${params.parent_ofall_dir}test/test12"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.NUMBER = 2



// Modules dependencie section

include { oneliner } from "${params.PIPES}fasta_oneliner" addParams(PUBLISH: "false")
include { counter } from "${params.PIPES}fasta_counter"




process split_number_computer {
	tag { "${fastas}" }
	container params.CONTAINER

	input:
	val tot_num_entries
	val desired_num_entries

	output:
	stdout emit: standardout

	script:
	"""
	#!/usr/bin/env python3
	
	#print("${tot_num_entries}")
	#print("${desired_num_entries}")
	list_entry = "${tot_num_entries}".split(',')
	tot_num_headers = 0
	list_of_tuples = []
	for elem in list_entry:
		#print(elem)
		ids_per_file = int((elem.split('=')[2].strip()).split(' ')[0])
		filename = ((elem.split('=')[3].strip()).split(',')[0]).split(']')[0]
		list_of_tuples.append((filename, ids_per_file))
		#print(ids_per_file)
		tot_num_headers += ids_per_file
	list_of_tuples.sort(key=lambda tup: tup[1], reverse=True)
	#print(tot_num_headers)
	#print(list_of_tuples)
	ratio = tot_num_headers // int("${desired_num_entries}")
	num_iterations_with_add_ratio = tot_num_headers - (ratio * int("${desired_num_entries}"))
	header_on_wich_change_pace = num_iterations_with_add_ratio * (ratio+1)
	#print(ratio, num_iterations_with_add_ratio, header_on_wich_change_pace)
	i = 0
	n = 0
	m = 0
	j = 0
	hits = 0
	counter_line = list_of_tuples[0][1]
	if num_iterations_with_add_ratio == 0:
		n = ratio
	else:
		n = ratio + 1
	#print("n =", n, "counter_line =", counter_line)
	for i in range(1,(tot_num_headers+1)):
		#print(i)
		if i == n:
			hits += 1
			j += 1
			#print("hits =", hits, "j =", j)
			#print("n =", n)
			if j < num_iterations_with_add_ratio:
				n += (ratio+1)
			else:
				n += ratio
			#print("n =", n)
		if i == counter_line:
			#print("i =", i)
			num_to_select = hits
			print(list_of_tuples[m][0], list_of_tuples[m][1], num_to_select, end=',')
			hits = 0
			if m < (len(list_of_tuples) - 1):
				m += 1
				counter_line += list_of_tuples[m][1]
	"""
}



process non_random_selector {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	tag { "${one_line_fastas}" }
	container params.CONTAINER

	input:
	path one_line_fastas
	val tot_num_entries_requested
	val desired_num_entries
	path py_script1

	output:
	path "*.${tot_num_entries_requested}sampledfasta", emit: sampled_fasta
	stdout emit: standardout

	script:
	out_name = "${one_line_fastas}".split('\\.')[0] + ".${tot_num_entries_requested}sampledfasta"
	list_tmp = "${desired_num_entries}".split(',')
	num_to_be_selected = ""
	for (item in list_tmp) {
		if ( item.contains("${one_line_fastas}") ) {
			num_to_be_selected = item
		}
	}
	"""
	echo ${one_line_fastas} 'n of headers ='`echo ${num_to_be_selected} | cut -d ' ' -f 2` 'selected ='`echo ${num_to_be_selected} | cut -d ' ' -f 3`
	./${py_script1} ${one_line_fastas} `echo ${num_to_be_selected} | cut -d ' ' -f 2` `echo ${num_to_be_selected} | cut -d ' ' -f 3` ${out_name}
	"""
}



workflow fasta_sampler {
	
	take:
	in_fastas
	number_requested

	main:
	non_random_py = params.SCRIPTS + "non_random_fasta_selecter.py"
	fasta_oneline = oneliner(in_fastas, "oneline")
	//fasta_oneline.view()
	//tmp = fasta_oneline.getClass()
	//println(tmp)
	tot_num_ids = counter(fasta_oneline)
	//tot_num_ids.view()
	split_number_computer(tot_num_ids.toList(), number_requested)
	ids_per_split = split_number_computer.out.standardout
	//ids_per_split.view()
	non_random_selector(fasta_oneline,  number_requested, ids_per_split.collect(), non_random_py)

	emit:
	sampledfasta = non_random_selector.out.sampled_fasta
	stout = non_random_selector.out.standardout
}



workflow {
	//count_out = fasta_counter(params.INPUT)
	fasta_sampler(params.INPUT, params.NUMBER)
	fasta_sampler.out.stout.view()
	fasta_sampler.out.sampledfasta.view()
	//println ( count_out.stout.getClass() )
}


//workflow.onComplete {
//	println ( workflow.success ? "\nDone! This is the output file bubba" )
//}
