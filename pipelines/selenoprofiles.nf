#!/usr/bin/env nextflow


// If the --help parameter is used in the command line, the pipeline will print the following help section:

if (params.help) {
	log.info 'This is the help section of this pipeline'
	log.info 'This pipeline performs a selenoprotein annotation of nucleotide sequences using selenoprofiles, which is a tool for predicting selenoprotein genes and their cognate SECIS elements. The pipeline can be run in two modes: (1) build a selenoprotein profile using multiple sequence alignments (MSA) as input, or (2) run the annotation using pre-built selenoprotein profiles.'
	log.info 'It can however be used to annotate other gene families with the correct profile.'
	log.info '\nFor more info about selenoprofiles4 functionig refer to the documentation page -> https://selenoprofiles4.readthedocs.io/en/latest/index.html\n'
	log.info 'The input parameters that manage the functioning of the pipeline are:'
	log.info '--CONTAINER		specifies the docker container image that contains the required dependencies and software versions for running the pipeline.'
	log.info '--INPUT		mandatory flag, specifies the input path to the nucleotide sequence file(s) in fasta format.' 
	log.info '			This are the file(s) that have to be predicted by selenoprofiles, a whole genome a chromosome or a gene region.'
	log.info '--INPUT_ALN		This parameter is optional, it switch to mode (2), it is the path to the multiple sequence alignment file(s) in fasta format.' 
	log.info '			Friom this file it builds custom selenoprotein profiles used during the DNA prediction step.'
	log.info '                      If more then one file is given as input each DNA file will be predicted with each custom profile.'
	log.info '--PROFILE		optional flag, it is mandatory if no ---INPUT_ALN has been given.'
	log.info '			It specifies the name of the pre-build profiles, for possible values refer to selenoprofiles4 documentation page.'
	log.info '--OUTPUT_DIR		optional flag, specifies the output directory for the selenoprotein annotation results.'
	log.info '			default value ${params.TEST_DIR}seleno_out/ , where the value of params.TEST_DIR is found in the nextflow.config file'
	log.info '--OUTPUT_FORMAT1	optional flag, specifies the output format for the predicted genes. default value   fasta  .'
	log.info '--OUTPUT_FORMAT2	optional flag, specifies the output format for the predicted genes. default value   gff  .'
	log.info '--PUBLISH		specifies whether to publish or not the output file and the custom profile, default true = publish = output'
	log.info '--SPECIES		optional flag, the name of the species analized, it can be whatever string value. It is used in internal filenames so be carefull.'
	log.info '			it is the value of the oprion -s in selenoprofiles.'
	log.info '\n####   WARNING  #####\n'
	log.info 'the pipeline always outputs 4 files + custom profile (3 files) if in mode(2). Out of the non-profile files 2 are always present:'
	log.info '.p2g and .ali (defaults in selenoprofiles4) on top of this the pipeline outputs two other files   .fasta and .gff  by default'
	log.info 'However this last 2 can be changed to any other type of supported output format. For the complete list refer to selenoprofiles4 documenattion page.'
        log.info '\n'
        exit 1
}



params.CONTAINER = "alessiovignoli3/tango-project:selenoprofiles" // selenoprofiles 4.4.3
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.INPUT_ALN = null
params.PROFILE = false
params.OUTPUT_DIR = "${params.TEST_DIR}seleno_out/"
params.OUTPUT_FORMAT1 = "fasta"
params.OUTPUT_FORMAT2 = "gff"
params.PUBLISH = true
params.SPECIES = "homo_sapiens"
params.BLAST_FILTER = 'default'
params.P2G_FILTER = 'default'
params.P2G_REFILTER = 'default'


process  generate_config_file {
	container params.CONTAINER
	
	output:
        path ".selenoprofiles_config.txt", emit: config_file, hidden: true
	
	script:
        """
        selenoprofiles -setup
        """
}



process apply_custom_filters {
        container params.CONTAINER
	stageInMode 'copy'

        input:
        path config_file
	val blast
	val filter
	val refilter

        output:
        path ".selenoprofiles_config.txt", emit: config_file, hidden: true
        stdout emit: standardout						// for debug

        script:
        """
	if [ ${blast} != 'default' ]; then
		sed -i "s/x.evalue < 1e-2  or x.sec_is_aligned()/${blast}/" .selenoprofiles_config.txt
	fi
	if [ ${filter} != 'default' ]; then
                sed -i "s/len(x.protein()) >60 or x.coverage()> 0.4/${filter}/" .selenoprofiles_config.txt
        fi
	if [ ${refilter} != 'default' ]; then
                sed -i "s/x.awsi_filter()/${refilter}/" .selenoprofiles_config.txt
        fi
        cat .selenoprofiles_config.txt
        """
}



process seleno_build_profile {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: true)
        tag { "${inaln}" }
        container params.CONTAINER

	input:
	path inaln
	path config_file

	output:
        tuple path("${out_profile_name}"), path("${out_profile_name}.config"), path("${out_profile_name}.profile_data"), emit: tupled_profile
        stdout emit: standardout

	script:
	out_profile_name = "${inaln.BaseName}" + "_profile.fa"
        """
        selenoprofiles build -i ${inaln} -o ${out_profile_name} -y
        """
}



process seleno_runner_custom {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true, saveAs:  { filename -> if(filename.endsWith(".ali")) "${out_names}.${filename.split('\\.')[-1]}"
										else "${out_names}.${filename.split('\\.')[-3]}.${filename.split('\\.')[-2]}.${filename.split('\\.')[-1]}" })
	tag { "${infasta}" }
	container params.CONTAINER
	//scratch true 

	input:
	each path(infasta)
	tuple path(profile), path(pro_config), path(pro_data)
	path config_file
	val species
	
	output:
	path "${profile.simpleName}*", emit: out_files, optional: true
	stdout emit: standardout

	script:
	out_names = "${infasta.simpleName}_" + "${profile.simpleName}"
	"""
 	selenoprofiles -o tmp -t ${infasta} -s ${species} -P ${profile} -output_${params.OUTPUT_FORMAT1} -output_${params.OUTPUT_FORMAT2}
	if [ "\$(ls -A tmp/${species}.${infasta.simpleName}/output/)"  ]; then
		mv tmp/${species}.${infasta.simpleName}/output/* .
	fi
	"""
}



process seleno_runner {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true, saveAs:  { filename -> if(filename.endsWith(".ali")) "${out_names}.${filename.split('\\.')[-1]}"
                                                                                else "${out_names}.${filename.split('\\.')[-3]}.${filename.split('\\.')[-2]}.${filename.split('\\.')[-1]}" })
	tag { "${infasta}" }
        container params.CONTAINER
	scratch true

        input:
        path infasta
	each profile
	path config_file
	val species

	output:
	tuple path("tmp/*/output/*.p2g"), path("tmp/*/output/*.ali"), emit: out_files, optional: true
        stdout emit: standardout

        script:
        """
        selenoprofiles -download -y
        selenoprofiles -o tmp -t ${infasta} -s ${species} -P ${profile} -output_${params.OUTPUT_FORMAT1}
	if [ "\$(ls -A tmp/${species}.${infasta.simpleName}/output/)"  ]; then
                mv tmp/${species}.${infasta.simpleName}/output/* .
        fi
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

	// Creating the config file only once at the beginning
	conf_file = generate_config_file()

	// Modifying the config file according to the filtering criterias if necessary
	config_with_filter = conf_file
	if (params.BLAST_FILTER != 'default' ||  params.P2G_FILTER != 'default' || params.P2G_REFILTER != 'default') {
		apply_custom_filters(conf_file, params.BLAST_FILTER, params.P2G_FILTER, params.P2G_REFILTER)
		config_with_filter = apply_custom_filters.out.config_file	
	}

	output_files = null
	if (pattern_to_aln != null) {
		in_aln = Channel.fromPath(pattern_to_aln)
		seleno_build_profile(in_aln, conf_file)
		seleno_runner_custom(in_fasta, seleno_build_profile.out.tupled_profile, config_with_filter, spcies_name)
		output_files = seleno_runner_custom.out.out_files
	} else {
		profiles = Channel.of("${profile_name}".split(","))
		seleno_runner(in_fasta, profiles, config_with_filter, spcies_name)
		output_files = seleno_runner.out.out_files
	}

	emit:
	outfile = output_files
	stout = apply_custom_filters.out.standardout		// for debug
}	


workflow {
	seleno_louncher(params.OUTPUT_DIR, params.INPUT, params.INPUT_ALN, params.PROFILE, params.SPECIES)
	seleno_louncher.out.outfile.view()
	seleno_louncher.out.stout.view()			// for debug
}

workflow.onComplete { println 'Done' }
