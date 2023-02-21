#!/usr/bin/env nextflow

// ### PLEASE CHECK THE COMMENTS AND COMMENTED LINES IN THE SCRIPTS or SHELLS BLOCKS ###




// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
      	log.info 'Here is a brief description of the pipeline itself  : '
        log.info '	this module takes as input a one id-keyword per line file and the fasta file to search such id and retrieve'
        log.info '	it takes care if multiple id files are given as input using the * as glob pattern same if many fasta are given'
        log.info '	but when both id files and fasta files are many a numeric keyword in the id nameand fasta names has to be given'
        log.info '	this means that each id file will be compared to one single fasta that has such keyword. For example:'
        log.info '	id_file_1.txt will be searched onto some_1.fasta the keyword has to be before the first point and after a _ character'
        log.info '	id_file_2.txt                       some_2.fasta'
        log.info '	or like this         some_otehr-id-file-glob_a.txt           other_fasta-glob_a.fasta.onlyfasta'
        log.info '\n'
        exit 1
}



params.OUTPUT_DIR = "${params.TEST_DIR}"
params.INHEADER = "${params.OUTPUT_DIR}test19_1.txt"
params.INFASTA =  "${params.OUTPUT_DIR}test18_1.fasta"
params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"	// slim buster



// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"



process retriever {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	tag { "${outname}" }
        container params.CONTAINER

	input:
	path inheader
	path fasta
	path py_script1

	output:
	path "${outname}", emit: outfasta
	stdout emit: standardout

	script:
	outname = ("${inheader}".split('\\.'))[0] + "-" + ("${fasta}".split('\\.'))[0] + ".fasta"
	"""
	python3 ${py_script1} ${inheader} ${fasta} ${outname} 2>tmp.err
        cat tmp.err
	"""

}



process pair_retriever {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	tag { "${outname}" }
	container params.CONTAINER

	input:
	tuple val(id), path(in_ids), path(in_fasta)
	path py_script1

	output:
	path "${outname}", emit: outfasta
	stdout emit: standardout

	script:
	outname = "${in_ids}".split('\\.')[0] + '-' + "${in_fasta}".split('\\.')[0] + '.fasta'
	"""
	./${py_script1} ${in_ids} ${in_fasta} ${outname} 2>tmp.err
	cat tmp.err
	"""
}




workflow fasta_retriever {

	take:
	pattern_to_idfile
        pattern_to_fastas

	main:
	in_idfiles = Channel.fromPath(pattern_to_idfile)
	in_fastas = Channel.fromPath(pattern_to_fastas)
	py_script_from = params.SCRIPTS + "from_keyword_to_fasta.py"
	//input_number_checcker(in_idfiles.count(), in_fastas.count())
	//check_line = path(input_number_checcker.out.standardout)
	splitted_in = ("${pattern_to_idfile}".split('\\*')).size()
	splitted_fa = ("${pattern_to_fastas}".split('\\*')).size()
	//println(splitted_in)
	//println(splitted_fa)
	//check_line.view()
	//println("${check_line}")
	if( splitted_in == 1 && splitted_fa > 1) {
		retriever(in_idfiles.collect(), in_fastas, py_script_from)
		tmp_stout = retriever.out.standardout
		tmp_outfile = retriever.out.outfasta
	} else 

	if( splitted_in > 1 && splitted_fa == 1) {
		retriever(in_idfiles, in_fastas.collect(), py_script_from)
		tmp_stout = retriever.out.standardout
		tmp_outfile = retriever.out.outfasta
	} else 

	if( splitted_in == 1 && splitted_fa == 1) {
		retriever(in_idfiles, in_fastas, py_script_from)
		tmp_stout = retriever.out.standardout
		tmp_outfile = retriever.out.outfasta
	} else {
		pairer(pattern_to_idfile, pattern_to_fastas)
		//pairer.out.right_pairs.view()
		pair_retriever(pairer.out.right_pairs, py_script_from)
		tmp_stout = pair_retriever.out.standardout
		tmp_outfile = pair_retriever.out.outfasta		
	}
		
	
	emit:
	//stout1 = check_line
	stout = tmp_stout
	retrieved = tmp_outfile
}



workflow {
	fasta_retriever(params.INHEADER, params.INFASTA)
	fasta_retriever.out.retrieved.view()
	//fasta_retriever.out.stout1.view()
	fasta_retriever.out.stout.view()
}



