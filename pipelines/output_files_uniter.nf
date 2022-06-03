#!/usr/bin/env nextflow

/*
####    ALERT    ####
#
# you might want to check collecFile() operator it does what is done here
#  there is an example in the script adjacent_annotation_computer.nf
*/

/*
*  This module is created for solving the issue with uniting into one single file the full content of 
*  many other files coming from the same process, it is usefull when a process is executed x times generating one file
*  each time and then the consequent pipeline needs onl one file that gathers what is inside them all.
*  
*/


nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      This module is created for solving the issue with uniting into one single file the full content of'
        log.info '      many other files coming from the same process, it is usefull when a process is executed x times generating one file'
        log.info '      each time and then the consequent pipeline needs only one file that gathers what is inside them all.'
        log.info '	it is made possible to laucìnch it alone for debug simplification'
        log.info '      It is necessary for the conteniriation to find the deepest common parent dir of all the files that have to be united'
        log.info '      so the first process is there only to solve this issue'
        log.info '      '
        log.info '      the --SWITCH_PUBLISH flag is set to false and it will not publish the result file, otherwise set it to whataerver'
        log.info '      you want'
        log.info '      the name of the output file is decided as <prefix>.<suffix>, this two variables are input parameters of the wokflow'
        log.info '              '
        log.info '              '
        exit 1
}

params.SWITCH_PUBLISH = false
params.INPUT1 = "${params.TEST_DIR}bubbabubba"
params.INPUT2 = "${params.TEST_DIR}bubbabubba"

process deepest_common_parentdir {
	//tag { "${outpath}" }
	container "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster

	input:
	val list_of_paths

	output:
	stdout emit: standardout

	script:
	"""
	#!/usr/bin/env python3
	
	#print(type('${list_of_paths}'))
	filepath_list = '${list_of_paths}'.split('[')[1].split(']')[0].split(', ')
	deepest_parentdir = filepath_list[0][::-1]
	for item in filepath_list:
		dirs_list = (item[::-1]).split('/')
		i = None
		for fol in dirs_list:
			i = deepest_parentdir.find(fol)
			if i != -1 :
				deepest_parentdir = deepest_parentdir[i:]
				break
		if i == -1:
			print('something went wrong not even root / has been found as common parent dir, check for sspaces and commas or other special charachters in the list of filpaths here below :')
			print('${list_of_paths}')
	print(deepest_parentdir[::-1] + '/')
	"""	
}

process single_channel_converger {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if ("${params.SWITCH_PUBLISH}" == "false") null
										else filename
										})
	tag { "${outname}" }
	container "ubuntu@sha256:86ac87f73641c920fb42cc9612d4fb57b5626b56ea2a19b894d0673fd5b4f2e9"
	
	input:
	val list_in_files
	path expedient
	val prefix
	val suffix
	
	output:
	path "${outname}", emit: converged
	stdout emit: standardout		// used in debug

	script:
	outname = "${prefix}.${suffix}"
	"""
	for i in `echo ${list_in_files}`; do cat `echo \$i | cut -d '[' -f 2 | cut -d ',' -f 1 | cut -d ']' -f 1` >> ${outname} ; done	
	echo ${outname}
	#echo ${list_in_files} ${expedient}
	"""
} 

workflow converger {
	
	take:
	input_channel
	prefix_name
	suffix_name

	main:
	deepest_common_parentdir(input_channel.toList())
	single_channel_converger(input_channel.collect(), deepest_common_parentdir.out.standardout, prefix_name, suffix_name)
	
	emit:
	stout = single_channel_converger.out.standardout
	converged_file = single_channel_converger.out.converged
}

workflow {
	in_channel = Channel.fromPath(params.INPUT1 + "*")
	converger(in_channel, 'test', 'bubba')
	converger.out.stout.view()
}
