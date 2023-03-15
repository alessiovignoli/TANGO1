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
params.OUTPUT_DIR = "${params.TEST_DIR}/seleno_out/"
params.PUBLISH = true
params.SPECIES = "homo_sapiens"


process seleno_runner {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true, saveAs: { filename -> if (params.PUBLISH == true) filename
										else null
										})	
	tag { "${infasta}" }
	container params.CONTAINER

	input:
	path output_folder
	path infasta
	val species
	

	output:
	//path "${output_folder}/*", emit: out_files
	stdout emit: standardout

	script:
	"""
	#sed -i '25i print("home config dir :", home_config_filename)' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py
	#sed -i '26i print("install dir ;", selenoprofiles_install_dir)' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py
	#sed -i '27i print("dependencies :", os.path.abspath(obonet.__file__),  os.path.abspath(easyterm.__file__),  os.path.abspath(sys.__file__))' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py
	#ls -alt /usr/local/selenoprofiles_workdir
	selenoprofiles -setup
	selenoprofiles -download -y
 	selenoprofiles -o tmp -t ${infasta} -s ${species} -P SelT,SelN,SPS
	"""
}


workflow seleno_louncher {
	
	take:
	patter_to_outdir
	pattern_to_fasta
	spcies_name

	main:
	out_dir = Channel.fromPath(patter_to_outdir, type: 'dir')
	in_fasta = Channel.fromPath(pattern_to_fasta)
	seleno_runner(out_dir, in_fasta, spcies_name)

	emit:
	outfile = 'bubba'
	stout = seleno_runner.out.standardout
}	


workflow {
	seleno_louncher(params.OUTPUT_DIR, params.INPUT, params.SPECIES)
	seleno_louncher.out.outfile.view()
	seleno_louncher.out.stout.view()
}


