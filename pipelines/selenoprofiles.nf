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
params.PROFILE = false
params.OUTPUT_DIR = "${params.TEST_DIR}seleno_out/"
params.PUBLISH = true
params.SPECIES = "homo_sapiens"



process seleno_build_profile {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: true)
        tag { "${inaln}" }
        container params.CONTAINER

	input:
	path inaln

	output:
        tuple path("${out_profile_name}"), path("${out_profile_name}.config"), path("${out_profile_name}.profile_data"), emit: tupled_profile
        stdout emit: standardout

	script:
	out_profile_name = "${inaln.BaseName}" + "_profile.fa"
        """
        selenoprofiles -setup
        selenoprofiles build -i ${inaln} -o ${out_profile_name} -y
        """

}



process seleno_runner_custom {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true)
	tag { "${infasta}" }
	container params.CONTAINER

	input:
	each path(infasta)
	tuple path(profile), path(pro_config), path(pro_data)
	val species
	
	output:
	tuple path("tmp/*/output/*.p2g"), path("tmp/*/output/*.ali"), emit: out_files, optional: true
	stdout emit: standardout

	script:
	"""
	selenoprofiles -setup
 	selenoprofiles -o tmp -t ${infasta} -s ${species} -P ${profile}
	"""
}



process seleno_runner {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true)
        tag { "${infasta}" }
        container params.CONTAINER

        input:
        path infasta
	each profile
	val species

	output:
	tuple path("tmp/*/output/*.p2g"), path("tmp/*/output/*.ali"), emit: out_files, optional: true
        stdout emit: standardout

        script:
        """
        selenoprofiles -setup
        selenoprofiles -download -y
        selenoprofiles -o tmp -t ${infasta} -s ${species} -P ${profile}
        """
}



workflow seleno_louncher {
	
	take:
	pattern_to_outdir
	pattern_to_fasta
	pattern_to_aln
	profile_name
	spcies_name

	main:
	in_fasta = Channel.fromPath(pattern_to_fasta)

	output_files = null

	if (pattern_to_aln != null) {
		in_aln = Channel.fromPath(pattern_to_aln)
		seleno_build_profile(in_aln)
		seleno_runner_custom(in_fasta, seleno_build_profile.out.tupled_profile, spcies_name)
		output_files = seleno_runner_custom.out.out_files
	} else {
		profiles = Channel.of("${profile_name}".split(","))
		seleno_runner(in_fasta, profiles, spcies_name)
		output_files = seleno_runner.out.out_files
	}

	emit:
	outfile = output_files
}	


workflow {
	seleno_louncher(params.OUTPUT_DIR, params.INPUT, params.INPUT_ALN, params.PROFILE, params.SPECIES)
	seleno_louncher.out.outfile.view()
}


