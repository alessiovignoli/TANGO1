#!/usr/bin/env nextflow

/*
*  This module is created for solving the issue with pairinf input files
*  it usefull in a situation where many couples of files that have to be paired but fromFilePairs is usable
*  for example there are three files in dir1 that have the names as file1.txt file2.txt file3.txt
*  and there is a correspondance with the three files in dir2 (or dir1 again) file1.out file2.out file3.txt
*  so file1.txt has to be processed only with file1.out and no the others 
*  same goes for file2 and 3 ecc..
*  this modules basically creates the correct input couple for other scripts that might need that
*  the file described so far are present in one channel and the second group in another
*  the script will take the two cHANNELS AS IMPUTS AND PAIR THE CORRESPONDING FILES INTO A SINGLE CHANNEL
*/


nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
	log.info 'This is the help section of this pipeline'
	log.info 'Here is a brief description of the pipeline itself  : '
	log.info '	This module is created for solving the issue with pairinf input files'
	log.info '	it usefull in a situation where many couplesof files that have to be paired but fromFilePairs is usable'
	log.info '      for example there are three files in dir1 that have the names as file1.txt file2.txt file3.txt'
	log.info '      and there is a correspondance with the three files in dir2 (or dir1 again) file1.out file2.out file3.txt'
	log.info '      so file1.txt has to be processed only with file1.out and no th4e others'
	log.info '      same goes for file2 and 3 ecc..'
	log.info '      this modules basically creates the correct input cuple for other scripts that might need that'
	log.info '      the file described so far are present in one channel and the second group in another'
	log.info '	the script will take the two cHANNELS AS IMPUTS AND PAIR THE CORRESPONDING FILES INTO A SINGLE CHANNEL'
	log.info '              '
	exit 1
}



params.INPUT1 = "${params.TEST_DIR}bubbabubba"
params.INPUT2 = "${params.TEST_DIR}bubbabubba"


workflow ch_pairer {
	
	take:
	input_files1
	input_files2

	main:
	input_files1.map{ it -> [it.simpleName, it]}.set{tupled_1}
	input_files2.map{ it -> [it.simpleName, it]}.set{tupled_2}
	//tupled_1.view()
	//tupled_2.view()
	tupled_1.combine(tupled_2, by:0).set{correct_pairs}
	//correct_pairs.view()
	correct_pairs.ifEmpty( { println('No matching scheme has been found between the input files. Remember it should be file1.some.txt file1.some.fasta,  <file1.some> is the matching key.')})
	
	emit:
	correct_pairs
}

workflow {
	ch1 = Channel.fromPath(params.INPUT1)
	ch2 = Channel.fromPath(params.INPUT2)
	ch_pairer(ch1, ch2)
}
