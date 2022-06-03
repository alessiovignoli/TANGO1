#!/usr/bin/env nextflow

/*
*  this module takes as input a txt file with the short ouptut prediction of phobius and
*  a corresponding fasta file, if in the txt file are present more than one sequences the programm expects
*  the same axact sequence in the same order as the txt in the  multifasta file
*  this module computes the amminoacid mass associated to whatever region of the prediciton is inputed
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
        log.info '      a corresponding fasta file option --INPUT_FATSA, if in the txt file are present more than one sequences the programm expects'
        log.info '      the same axact sequence in the same order as the txt in the fasta file, multifasta file'
        log.info '      this module computes the amminoacid mass associated to whatever region of the prediciton is inputed'
        log.info '      this region can be selected simply giving in the command line the  --KEYWORD argument'
        log.info '\n'
        exit 1
}



params.CONTAINER = "alessiovignoli3/tango-project:python_field_retr@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"
params.INPUT_TXT = "${params.TEST_DIR}bubbabubba"
params.INPUT_FATSA = "${params.TEST_DIR}bubbabubba"
params.KEYWORD = false
params.MAX_ITER = false
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.OUTPUT_NAME = "bubba"



// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"
include { average_unifier } from "${params.PIPES}average_plp_short_pred" addParams(SUFFIX: "aacomp") 


process based_on_short_aamass_fasta {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if (params.MAX_ITER == false) filename
	//									else  null
	//									})
	tag { "${fasta}" }
	container params.CONTAINER

	input:
	tuple val(id), path(txt), path(fasta)
	path pyscript1
	val field_id

	output:
	//path "${output_name}", emit: aamass_fasta
	stdout emit: standardout

	script:
	prefix_fasta = "${fasta}".split('\\.')[0]
	output_name = prefix_fasta + ".aamass"
	"""
	./${pyscript1} ${txt} ${fasta} ${field_id} ${params.MAX_ITER} #2>tmp.err 1>${output_name}
	#cat tmp.err
	"""
}



workflow short_pred_amminoacid_mass_fasta {
	
	take:
	input_txt
	input_fasta
	field_keyword

	main:
	fasta_mass_computer_py = params.SCRIPTS + "fasta_short_pred_hydrophobicity_computer.py"
	pairer(input_txt, input_fasta)
	//pairer.out.right_pairs.view()
	based_on_short_aamass_fasta(pairer.out.right_pairs, fasta_mass_computer_py, field_keyword)
	//aamassfasta = based_on_short_aacomp_fasta.out.aamass_fasta
	finalaacompfasta = "there is nothing here, the --MAX_ITER has not been given a values, hence no need for this block to be executed this is not an error"
	if ( params.MAX_ITER != false) {
		average_unifier(aamassfasta.collect(), fasta_aamass_computer_py)
		finalmassfasta = average_unifier.out.final_average
	}

	emit:
	//aamassfasta
	//finalaacompfasta
	stout1 = based_on_short_aamass_fasta.out.standardout
	//stout2 = average_unifier.out.standardout
}



workflow {
	short_pred_amminoacid_mass_fasta(params.INPUT_TXT, params.INPUT_FATSA, params.KEYWORD)
	short_pred_amminoacid_mass_fasta.out.stout1.view()
	//short_pred_amminoacid_mass_fasta.out.stout2.view()
	//short_pred_amminoacid_mass_fasta.out.aacompfasta.view()
	//short_pred_amminoacid_mass_fasta.out.finalaacompfasta.view()	
}
