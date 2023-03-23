#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info ''
        log.info ''
        log.info ''
        log.info '\n'
        exit 1
}



params.CONTAINER = "alessiovignoli3/tango-project:selenoprofiles" // selenoprofiles 4.4.3
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.INPUT_ALN = null
params.OUT_PROFILE = false
params.OUTPUT_DIR = "${params.TEST_DIR}seleno_out/"
params.PUBLISH = true
params.SPECIES = "homo_sapiens"



process seleno_build_profile {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true)
        tag { "${inaln}" }
        container params.CONTAINER

	input:
	path inaln
	val out_profile_name

	output:
        //path "tmp/*", emit: out_dir_files
        stdout emit: standardout

	script:
	if (out_profile_name == false) {
		out_profile_name = "${inaln.BaseName}" + "_profile.fa"
	}
        """
        selenoprofiles -setup
        #selenoprofiles -download -y
        selenoprofiles build -i ${inaln} -o ${out_profile_name} -y
        """

}



process seleno_runner {
	//publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true)
	tag { "${infasta}" }
	container params.CONTAINER

	input:
	path infasta
	val species
	

	output:
	//path "tmp/**", emit: out_dir_files
	stdout emit: standardout

	script:
	"""
	selenoprofiles -setup
	selenoprofiles -download -y
 	selenoprofiles -o tmp -t ${infasta} -s ${species} -P SelT,SelN,SPS
	"""
}


workflow seleno_louncher {
	
	take:
	pattern_to_outdir
	pattern_to_fasta
	pattern_to_aln
	outprofile_name
	spcies_name

	main:
	//out_dir = Channel.fromPath(pattern_to_outdir, type: 'dir')
	in_fasta = Channel.fromPath(pattern_to_fasta)
	
	if (pattern_to_aln != null) {
		in_aln = Channel.fromPath(pattern_to_aln)
		seleno_build_profile(in_aln, outprofile_name)
	}

	//seleno_runner(in_fasta, spcies_name)

	emit:
	outfile = 'bubba' //patter_to_outdir
	stout = seleno_build_profile.out.standardout //seleno_runner.out.standardout
}	


workflow {
	seleno_louncher(params.OUTPUT_DIR, params.INPUT, params.INPUT_ALN, params.OUT_PROFILE, params.SPECIES)
	seleno_louncher.out.outfile.view()
	seleno_louncher.out.stout.view()
}


