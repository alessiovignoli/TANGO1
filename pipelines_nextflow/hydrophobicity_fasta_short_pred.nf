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
*  this module takes as input a short prediction phobius output file (simple redirection of stdout) 
*  that has the following structure (space separated)
*	UniRef50_A0A5J4NZG1             2  Y n15-25c30/31o1387-n-1403i1415-n-1437o
*	UniRef50_W7A840                 1  0 o1786-n-1803i
*  where each line is the prediction for a specifi protein
*  it also takes as input the fasta file that generated such prediction (it needs to have the sequence on one line)
*  check the fasta_oneliner.nf script if you need
*  So this pipeline based on a segment obtained from the predfile computes the Hydrophobicity associated to such segment
*  for this reason the order of the proteins in the predfile must be the same for the fasta one
*  basically first protein in pred must also be the first protein in fasta and so on
*/



nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
      	log.info 'Here is a brief description of the pipeline itself  : '
        log.info '	this module takes as input a short prediction phobius output file --INPUT_TXT (simple redirection of stdout)'
        log.info '	that has the following structure (space separated)'
        log.info '		UniRef50_A0A5J4NZG1             2  Y n15-25c30/31o1387-n-1403i1415-n-1437o'
        log.info '		UniRef50_W7A840                 1  0 o1786-n-1803i'
        log.info '	where each line is the prediction for a specifi protein'
        log.info '	it also takes as input the fasta file that generated such prediction --INPUT_FASTA '
        log.info '	(it needs to have the sequence on one line) check the fasta_oneliner.nf script if you need'
        log.info '	So this pipeline based on a segment obtained from the predfile computes the Hydrophobicity associated to such segment'
        log.info '	for this reason the order of the proteins in the predfile must be the same for the fasta one'
        log.info '	basically first protein in pred must also be the first protein in fasta and so on'
        log.info '	both the previous flags are mandatory and ihf the input is a pattern to many files the script assumes that for each txt'
        log.info '	there is a corresponding fasta file, remember that to be recognized as pattern do like --INPUT_FASTA "test*.fasta"'
        log.info '	it is mandatory to specify the region of the prediciton on which do the computation, --KEYWORD argument '
        log.info '	allowed values are the allowed keywords are: [c, i, o, -, n, s, l]  c = signal peptide, i = inside membrane(cytoplasm),'
        log.info '      o = outside membrane, - = helix (in phobius originalmodel), (only in phobius-M7or later) => -n- = normal-helix'
        log.info '      -s- = special-helix and -l- = loop-inramembrane'
        log.info '      they have to be given like this:'
        log.info '              --KEYWORD s,n,c    or  i,l,\\-'
        log.info '      optionally the name for the hydrophobicity scale can be specified using --HYDRO_SCALE'
        log.info '	the supported arguments are kyte for Kyte-Doolittle (default), GES for GES scale, UHS for Unified Hydrophob Scale'
        log.info '	'
        log.info '\n'
        exit 1
}


params.CONTAINER = "alessiovignoli3/tango-project:python_field_retr@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.INPUT_TXT = "${params.OUTPUT_DIR}test19_1.txt"
params.INPUT_FASTA =  "${params.OUTPUT_DIR}test18_1.fasta"
params.KEYWORD = false
params.MAX_ITER = false
params.HYDRO_SCALE = "kyte"


// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"
include { average_unifier } from "${params.PIPES}average_plp_short_pred" addParams(INPUT_PLP: "${params.INPUT_TXT}", SUFFIX: "${params.HYDRO_SCALE}_avghydro")



process pred_avg_hydro_fasta {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false , saveAs: { filename -> if (params.MAX_ITER == false) filename
										else  null
										})
	tag { "${outname}" }
	container params.CONTAINER

	input:
	tuple val(id), path(in_txt), path(in_fasta)
	val field_id
	path py_script1

	output:
	path "${outname}", emit: avg_hydro
	stdout emit: standardout

	script:
	outname = "${in_txt}".split('\\.')[0] + "_${field_id}.${params.HYDRO_SCALE}_averagehydro"
	"""
	./${py_script1} ${in_txt} ${in_fasta} ${field_id} ${params.HYDRO_SCALE} ${params.MAX_ITER} 2>tmp.err 1>${outname}
	cat tmp.err 
	"""
}




workflow average_hydrophobicity {

	take:
	pattern_to_txt
        pattern_to_fastas
	field_keyword

	main:
	pairer(pattern_to_txt, pattern_to_fastas)
	//pairer.out.right_pairs.view()
	hydrophobicity_computer_py = params.SCRIPTS + "fasta_short_pred_hydrophobicity_computer.py"
	pred_avg_hydro_fasta(pairer.out.right_pairs, field_keyword, hydrophobicity_computer_py)
	avghydro = pred_avg_hydro_fasta.out.avg_hydro
	finalavghydro = "there is nothing here, the --MAX_ITER has not been given a values, hence no need for this block to be executed this is not an error"
	if ( params.MAX_ITER != false) {
		average_unifier(avghydro.collect(), hydrophobicity_computer_py, field_keyword)
		finalavghydro = average_unifier.out.final_average
	}
	
	emit:
	avghydro
	finalavghydro
	stout1 = pred_avg_hydro_fasta.out.standardout
}



workflow {
	average_hydrophobicity(params.INPUT_TXT, params.INPUT_FASTA, params.KEYWORD)
	//average_hydrophobicity.out.avghydro.view()
	average_hydrophobicity.out.finalavghydro.view()
	average_hydrophobicity.out.stout1.view()
}


