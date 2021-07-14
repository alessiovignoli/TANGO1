#!/usr/bin/env nextflow

/*
*  this module takes as input a plp file with otput of phobius with the following characteristics:
*  it can be a single sequence plp file or a multisequence multiplp file where every sequence prediction plp is separated by two line like this :
*       # sp|Q8BI84|TGO1_MOUSE
*       #pos    aa      i       o       M       S       m       l
*  where the first comment line contains the id of the sequence that generated such prediction, and the second line describes the content
*  of the columns
*  the output file will be similar to the plp or multiplp file, it will have this structure:
*  	protein_id	position_of_aa_in_sequence 	amminoacid_code		plp_cytopplasmic  plp_non_cytoplasmic  plp_normal_membrane  plp_signal_peptide  plp_special_helix  plp_loop_intramembrane
*	sp|Q8BI84|TGO1_MOUSE		1		00000000010000000000	0.00128		0.99872		0.00000		0	0.00000		0.00000
*  where the first line that describes the columns is only at the first line.
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      this module takes as input a plp file option --INPUT_PLP otput of phobius with the following characteristics:'
        log.info '      it can be a single sequence plp file or a multisequence multiplp file where every sequence prediction plp is separated by two line like this :'
        log.info '      	# sp|Q8BI84|TGO1_MOUSE'
        log.info '      	#pos    aa      i       o       M       S       m       l'
        log.info '      where the first comment line contains the id of the sequence that generated such prediction, and the second line describes the content'
        log.info '      of the columns'
        log.info '      the output file will be similar to the plp or multiplp file, it will have this structure:'
        log.info '      protein_id      position_of_aa_in_sequence      amminoacid_code         plp_cytopplasmic  plp_non_cytoplasmic  plp_normal_membrane  plp_signal_peptide  plp_special_helix  plp_loop_intramembrane'
        log.info '      sp|Q8BI84|TGO1_MOUSE            1               00000000010000000000    0.00128         0.99872         0.00000         0       0.00000         0.00000'
        log.info '      where the first line that describes the columns is only at the first line.'
        log.info '\n'
        exit 1
}


params.CONTAINER = "alessiovignoli3/tango-project:python_field_retr@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b"
params.INPUT_PLP = "${params.TEST_DIR}bubbabubba"
params.MAX_ITER = false
params.KEYWORD = false
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.OUTPUT_NAME = "bubba"


process aa_from_plp_hard_encoder {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	tag { "${average_outfile}" }
	container params.CONTAINER

	input:
	path plp
	path pyscript1

	output:
	path "${output_name}", emit: encoded_plp
	stdout emit: standardout

	script:
	output_name = "${plp}".split('\\.')[0] + "-${params.MAX_ITER}.endedplp"
	if ("${params.MAX_ITER}" == false ) {
		output_name = "${plp}".split('\\.')[0] + ".endedplp"
	}
	"""
	./${pyscript1} ${plp} ${params.KEYWORD} ${params.MAX_ITER} 2>tmp.err 1>${output_name}
	cat tmp.err
	"""
}



workflow aa_hard_encoder {

	take:
	input_plp
	
	main:
	in_plp = Channel.fromPath(input_plp)
	aa_encoder_py = params.SCRIPTS  + "plp_aa_encoder.py"
	aa_from_plp_hard_encoder(in_plp, aa_encoder_py)

	emit:
	out_plp = aa_from_plp_hard_encoder.out.encoded_plp
	stout1 = aa_from_plp_hard_encoder.out.standardout
}


workflow {
	aa_hard_encoder(params.INPUT_PLP)
	aa_hard_encoder.out.stout1.view()
	aa_hard_encoder.out.out_plp.view()
}
