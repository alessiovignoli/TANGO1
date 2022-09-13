#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline execute the uagustus ppx function, basically a prot-MSA refined gene prediction"
        log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info '--IN		file or filepath as glob pattern to genome files containing portion of DNA sequence in fasta format'
        log.info "		when many files are given, each one is treated as separate genome/specie"
        log.info "--REF_FA	give this flag if the sequences comprising the MSA used for enhanche the genome predicitons have to be"
	log.info "		alligned still. This flags accepts only one fasta file, so if sequences nedded for alignment are in differnt"
        log.info "		files is uop the user to put them in one. The produced alignment is going to be published in the OUT_DIR directory"
        log.info "--REF_ALN	pass this flag if the sequences are already alligned, the format accepted by augustus are: FASTA or CLUSTAL."
	log.info "              This flag accepts only one file."
	log.info "		"
	log.info "the above two flags are mutually exclusive, only one of them can be given."
        log.info "##########   WARNING    ###########"
	log.info "if more than one file is passed in one of the above --REF flags, for example as a glob path"
	log.info "only the first file is going to be used"
        log.info ""
        log.info "--PREP_ALN	optional flag, default false no neeed for   prepareAlign   module of Augustus."
        log.info "		basically if the input MSA or the one generated is a bit too gappy or messy or big"
	log.info "		then is better to give this flag and use such module specifically made for circumventing the problem."
        log.info "		More info can be found at    https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md"
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

params.CONTAINER = "alessiovignoli3/tango-project@sha256:9a351679d2f41c54b2baabb44feb415c870344406f11bde627854d97f810aaf9" // augustus 3.4.0
params.IN = "${params.TEST_DIR}bubbabu"
params.REF_FA = false
params.REF_ALN = false
params.OUT_DIR = "${params.TEST_DIR}augustus_results/"
params.PREP_ALN = false

// Modules dependencie section

include { oneliner_ch } from "${params.PIPES}fasta_oneliner" addParams(PUBLISH: "false", CONTAINER: "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43")	// the same used in this module
include { align_generation } from "${params.PIPES}aln_generation_colouring" addParams(OUTPUT_DIR: params.OUT_DIR) 


process  aug_prep_aln {
	container params.CONTAINER
	tag { "${aln}" }

	input:
	path aln

	output:
	//path "*", emit: tmp_files
	stdout emit: standardout	// for debug

	script:
	outname = "${aln}".split('\\.')[0] + '.preppedaln'
	"""
	prepareAlign < ${aln} > ${outname} 
	"""
}





workflow  augustus_ppx {

	take:
	pattern_to_genomes
	pattern_to_fasta
	pattern_to_msa

	main:

	// error section of missing or wrong inputs

	if ( !pattern_to_fasta && !pattern_to_msa ) {
		log.info "ERROR: nor --REF_FA nor --REF_ALN has been given, one of the two is mandatory, for more details type --help"
	} 

	
	// alignment preparation section

	in_ref = null
	if ( pattern_to_fasta ) {
		fastas = channel.fromPath(pattern_to_fasta)
		oneliner_ch(fastas.first(), 'onelinefa')
		align_generation(oneliner_ch.out.onelinefasta)
		in_ref = align_generation.out.aln_file
	} else {
		in_ref = channel.fromPath(pattern_to_msa)
	}
	// change from clustal to fasta allign, implement it in other file
	if ( params.PREP_ALN ) {
		aug_prep_aln(in_ref)
	}


	emit:
	stout = aug_prep_aln.out.standardout //preppedaln()
}

workflow {
	augustus_ppx(params.IN, params.REF_FA, params.REF_ALN)
	augustus_ppx.out.stout.view()
}

