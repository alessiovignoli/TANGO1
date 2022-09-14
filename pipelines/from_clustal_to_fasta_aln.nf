#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline converts a clustal Multiple Sequence allignment (MSA) into a fasta format alignment unsing tcoffee"
        log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info "--ALN		The clustal MSA file or filepath, it can be a glob path to multiple MSAs"
        log.info ""
        log.info "--SUFFIX	optional field, the string to attach at each MSA name file, default   -aln.fasta   "
        log.info "		the string given will be added at the end of the name in place of the extension."
        log.info "		if the filene is   something.whatever.fa   the output file will be called by default"
        log.info "		something-aln.fasta"
	log.info "--OUTPUT_DIR	optional flag, the directory where all output files are stored in, default params.TEST_DIR "
        log.info "		variable found in nextflow.config, it should be in the test/ dir"
        log.info "--PUBLISh	optional falg, default true, means put the output files in the directory specified by --OUTPUT_DIR."
        log.info "		usually set to false in case of module import in other pipelines."
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

params.CONTAINER = "cbcrg/tcoffee@sha256:36c526915a898d5c15ede89bbc3854c0a66cef22db86285c244b89cad40fb855"  // tcoffee 13.45.47.aba98c5 
params.ALN = false
params.SUFFIX = '-aln.fasta'
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.PUBLISH = true



process clustal_converter  {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: false, saveAs: { filename -> if (params.PUBLISH == true) filename
                                                                                else null
                                                                                })
	container params.CONTAINER
	tag { "${aln}" }

	input:
	path aln
	val extension 

	output:
	path "*${extension}", emit: fasta_aln
	stdout emit: standardout		// for debug

	script:
	outname = "${aln}".split('\\.')[0] + extension
	"""
	t_coffee -other_pg seq_reformat -in ${aln} -output fasta_aln > ${outname}
	"""
}



workflow clustal_to_fasta {

	take:
	pattern_to_aln
	suffix

	main:

	// error section

	if ( !pattern_to_aln ) {
		log.info "ERROR: the --ALN argument has not been given, type --help for more details"	
	}

	alns = channel.fromPath(pattern_to_aln)
	clustal_converter(alns, suffix)

	emit:
	aln_fasta = clustal_converter.out.fasta_aln
}

workflow {
	clustal_to_fasta(params.ALN, params.SUFFIX)
	clustal_to_fasta.out.aln_fasta.view()
}

