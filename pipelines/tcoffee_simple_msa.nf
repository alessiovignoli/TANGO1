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



params.CONTAINER =  "cbcrg/tcoffee@sha256:36c526915a898d5c15ede89bbc3854c0a66cef22db86285c244b89cad40fb855" // t_coffee Version_13.45.47.aba98c5
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTNAME = false
params.FORMAT = "clustalw,fasta_aln"
params.OUTPUT_DIR = "${params.TEST_DIR}seleno_out/"
params.PUBLISH = true




process simple_msa {
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: true)
        tag { "${infasta}" }
	container params.CONTAINER

	input:
        path infasta
	val outname
	

	output:
	path "*", emit: alns
	stdout emit: standardout

	script:
	out_file = outname
	if (outname==false) {
		"""
		t_coffee -in ${infasta}  -output ${params.FORMAT}
		rm -f *.dnd
		"""
	} else {
        	"""
		t_coffee -in ${infasta} -outfile ${out_file} -output ${params.FORMAT}
		rm -f *.dnd
		"""
	}
}





workflow tcoffee_simple_louncher {
	
	take:
	pattern_to_fasta
	pattern_outname

	main:
	in_fasta = Channel.fromPath(pattern_to_fasta)

	simple_msa(in_fasta, pattern_outname)
	

	emit:
	errors = simple_msa.out.standardout
	outfile = simple_msa.out.alns
}	


workflow {
	tcoffee_simple_louncher(params.INPUT, params.OUTNAME)
	tcoffee_simple_louncher.out.errors.subscribe onError: {println it }
	tcoffee_simple_louncher.out.outfile.subscribe onNext: { println it }, onComplete: { println 'Done' }
}


