#!/usr/bin/env nextflow

// ### PLEASE CHECK THE COMMENTS AND COMMENTED LINES IN THE SCRIPTS or SHELLS BLOCKS ###


// tm = tranmembrane
// om = Original model
// M3 = mark three third model created of the phobius model
// ns = negatie set
// ps = positive set TANGO1 proteins not used for the training of the model
// fs = frequency setthe TANGO1 proteins used for the assestment of the frequency of aa this prot are excluded from the positive set
// pred = prediction


nextflow.enable.dsl=2

params.help = false


/*
*  This module takes as input fasta files either single fasta or multifasta files and lauch phobius in a parallel manner
*  a it laounches a process for each input fasta file 
*  it has different options and is associated with some config profiles such as mark7_models.config, original_model.config, test.config
*  take a look at them in the conf dir (usually underthe test dir) the nextflow.config that launches such profiles is in the test dir
*/


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '      This module takes as input fasta files either single fasta or multifasta files and lauch phobius in a parallel manner'
        log.info '	--INPUT   option'
        log.info '      it laounches a process for each input fasta file'
        log.info '      it has different options and is associated with some config profiles such as mark7_models.config, original_model.config, test.config'
	log.info '	take a look at them in the conf dir (usually underthe test dir) the nextflow.config that launches such profiles is in the test dir'
        log.info '      '
        log.info '      The other two flags that can be passed to this pipeline are switches that rergulate the behaviour of the pipeline itself'
        log.info '      skipping some process or choosing between different workflows. The ffirst one is --ONE_LINE that takes as values the boolean true and false'
        log.info '      default false, and it means that the fasta files do not have the sequence on one line, the pipeline will put them on one line'
        log.info '      if the file are already with the sequence on one line pass an argumrnt to --ONE_LINE'
        log.info '      --PLP will prompt the pipeline to compute a plp file for each fasta sequence in the input'
        log.info '      default false it means no plp will be made'
        log.info '      '
        log.info '\n'
        exit 1
}



// params used in this script

params.infasta_paral = "bubba"
//params.CONTAINER = "alessiovignoli3/tango-project@${params.sha_code}" "alessiovignoli3/tango-project:phobius_image@${params.sha_code}"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.INPUT = "${params.TEST_DIR}${params.infasta_paral}"
params.ONE_LINE = false
params.PLP = false


// Modules dependencies section

include { oneliner } from "${params.PIPES}fasta_oneliner" addParams(OUTPUT_DIR: "/dev/null")




process phobius_short_parallelization {
	label 'short'
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)

	input:
	path fasta

	output:
	path "*txt", emit: predfile

	script:
	suffix = "${params.model_version}".split('phobius')[1]
	prefix = "${fasta}".split('fasta')[0]
	"""
	export PATH="/home/phobius101_linux/tmp/tmpbyom7j/phobius:$PATH"
	#tail -n 100 /home/phobius101_linux/tmp/tmpbyom7j/phobius/phobius.pl
	limit=`cat $fasta | wc -l`
	i=2; while [ \$i -le \$limit ]; do 
		SEQUENCE=`sed -n "\$i"p ${fasta}`; if echo "\$SEQUENCE" | grep -vq  "O"; then
    				sed -n  `expr \$i - 1`,"\$i"p ${fasta} | ${params.exec_version} -short
			fi;i=\$((\$i + 2)); done 1>"phobius_output${suffix}_${prefix}txt"
	"""
}


process phobius_short_and_plp_parallelization {
	label 'short_and_plp'
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	//container "alessiovignoli3/tango-project:splp_phobius_image@sha256:f098f1511f37d461f3610884c8197814231eeaec3d09286a941292fa5f289dd5"

	input:
	path fasta

	output:
	path "*txt", emit: predfile
	path "*plp", emit: plpfile
	
	script:
	suffix = "${params.model_version}".split('phobius')[1]
	prefix = "${fasta}".split('fasta')[0]
	"""
	export PATH="/home/phobius101_linux/tmp/tmpbyom7j/phobius:$PATH"
	#tail -n 100 /home/phobius101_linux/tmp/tmpbyom7j/phobius/phobius.pl
	#cat /home/phobius101_linux/tmp/tmpbyom7j/phobius/phobius_model_versions/${params.model_version} > /home/phobius101_linux/tmp/tmpbyom7j/phobius/phobius.model
	limit=`cat $fasta | wc -l`
	i=2; while [ \$i -le \$limit ]; do
		SEQUENCE=`sed -n "\$i"p ${fasta}`; if echo "\$SEQUENCE" | grep -vq  "O"; then
				sed -n  `expr \$i - 1`,"\$i"p ${fasta} | ${params.exec_version} -short -plp "${prefix}tmpplp"
				cat "${prefix}tmpplp" >> "phobius_output${suffix}_${prefix}plp"
			fi;i=\$((\$i + 2)); done 1>"phobius_output${suffix}_${prefix}txt"
	rm "${prefix}tmpplp"
        """
	
}

workflow phobius_short {

	take:
	input_fasta

	main:
	if (params.ONE_LINE != false) {
                in_fasta = Channel.fromPath(input_fasta)
                phobius_short_parallelization(in_fasta)
	} else {
		suffix = "oneline"
		fasta_oneline = oneliner(input_fasta, suffix)
		//fasta_oneline.view()
		phobius_short_parallelization(fasta_oneline)
	}
	
	emit:
	pred_file = phobius_short_parallelization.out.predfile
}


workflow phobius_short_plp {

	take:
	input_fasta

	main:
	if (params.ONE_LINE != false) {
		in_fasta = Channel.fromPath(input_fasta)
		phobius_short_and_plp_parallelization(in_fasta)
	} else {
		suffix = "oneline"
        	fasta_oneline = oneliner(input_fasta, suffix)
		phobius_short_and_plp_parallelization(fasta_oneline)
	}

	emit:
	pred_file = phobius_short_and_plp_parallelization.out.predfile
	plp_file = phobius_short_and_plp_parallelization.out.plpfile
}


workflow {
	if (params.PLP != false) {
		phobius_short_plp(params.INPUT)
        	phobius_short_plp.out.pred_file.view()
	        phobius_short_plp.out.plp_file.view()
	} else {
		phobius_short(params.INPUT)
		phobius_short.out.pred_file.view()
	}
}


