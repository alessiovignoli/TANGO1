#!/usr/bin/env nextflow




// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
	log.info "This pipeline is used to downlad from the HTTP server the ensmble info requested"
	log.info ""
	log.info ""
	log.info "Here is the list of flags accepted by the pipeline:"
	log.info "--INPUT_IDS		the file or files containig the query ENSEMBL IDs on one column, the file has to have all Ids on same column"
	log.info "			and the columns have to be separated always by the same field separator, specifiable, read below"
	log.info "--TYPE_INPUT		optional flag that specifies to the pipeline if input IDs are gene Id, transcript Id or protein ID."
        log.info "			default G for gene, T for transcript and P for protein, this are the values accepted"
        log.info "--TYPE_OUTPUT		optinal flag, specifies what ahs to be recovered from ENSEMBL, possible values: genomic, cds, cdna, protein"
        log.info "			default genomic (DNA not spliced, raw), for more detail take a look at the link at the bottom"
        log.info "#######     WARNING     #######"
        log.info "not all the combination of the above flags have an output or are already implemented, read below"
        log.info ""
	log.info "--FIELD_SEP		optional flag, tells the script which is the field/column separator used in the input file/s"
	log.info "			default \\t <tab>, to give something different pass --FIELD_SEP ';\\t' , in this example the columns are"
	log.info "			separated by a comma followed by a <tab>, quotation are necessary in this case. The extraction of the columns"
	log.info "			is done through awk program, specifically FS= variable, look at awk documentation for special charachtes"
	log.info "			for example field are separated by two or more spaces can be passed like this ' {2,}' or '  +' "
	log.info "--COL			optional flag, specifies which is the column to select, default 1 meaning first column"
	log.info "			for last column, usefull when the number of columns per line is variable, look at awk NF, give it as"
        log.info "			--column_id_pdb 'NF' or 'NF-1' for second to last."
        log.info "--MASK_INTRON		optional flag to mask intron when input and output are genomic sequences, default 0 no mask, give 1 for masking"
        log.info "			intron will be writen in lower case in that scenario"
        log.info "--OUTPUT_DIR		optional flag, the directory where all output files are stored in, default params.TEST_DIR variable found in"
        log.info "			nextflow.config, it should be the test/ dir"
        log.info ""
	log.info ""
	log.info "#######     WARNING     #######"
	log.info "is up to the user making sure the inoput IDs are correct and is also preferrable they have to be of the same type, example:"
        log.info "all Ensemb gene IDs ENSG0... or all protein IDs ENSP0... ecc.., this is due to the fact that some flags are only appliable to"
	log.info "gene Ids for example and will create errors on other type of IDs."
	log.info "Also keep in mind that ID like ENSAMET00000029331.1 will not produce any output, just remove the .1"
        log.info ""
        log.info "Also if thousands of files have to be downloaded better split the ID file (--query_pdb) into many smaller ones"
	log.info "so that the process is parallelized. The structure of the 'smaller' files have to be consistemt. If many ID input file have to be passed"
        log.info 'do like: --query_pdb "id_file*.tab" where there are id_file1.tab, id_file2.tab ecc..'
	log.info ''
        log.info ""
        log.info ""
        log.info "---- OUTPUT EXPLANATION ----"
        log.info "Here is a brief description of the combination of the two optional flags TYPE_INPUT and TYPE_OUTPUT and the results they generate"
        log.info "TYPE_INPUT		TYPE_OUTPUT"
        log.info "G			genomic		input gene IDs, output is the whole unspliced genomic region, it can be masked for introns"
        log.info "G			protein		as above, out are all the protein sequences (aa) associated whith such gene in xml format"
        log.info "T			cds		in transcript ID, out  the spliced transcript sequence without UTR in fasta format"
        log.info "T			cdna		as above, out spliced transcript sequence with UTR in fasta format"
        log.info "P			protein		input protein ID, out protein sequence in json format"
        log.info ""
        log.info ""
	log.info 'for nore details look at ->  https://rest.ensembl.org/documentation/info/sequence_id'
	log.info ''
	exit 1
}


params.CONTAINER = "alessiovignoli3/tango-project@sha256:57013bf372b519245608c95fd60a38f9e5d65775aaf18c2d711031818c1a145e"     // bash5.0.17 with awk and wget
params.INPUT_IDS = "${params.TEST_DIR}bubbabubba.ids"
params.TYPE_INPUT = "G"
params.TYPE_OUTPUT = "genomic"
params.FIELD_SEP = '\t'
params.COL = 1
params.MASK_INTRON = 0
params.OUTPUT_DIR = "${params.TEST_DIR}"


process get_ensembl_emtry {
	beforeScript "mkdir -p ${params.OUTPUT_DIR}"
	publishDir(params.OUTPUT_DIR, mode: 'move', overwrite: false)
        tag { "${ID_file}" }
        label 'pdb'
	container params.CONTAINER

        input:
        path ID_file
	val type_of_input
        val type_of_output

        output:
        path "*.fasta*", emit: fasta
        stdout emit: standardout                //for debug

        script:
        if (type_of_input=='G' && type_of_output=='genomic') {
                """
                for i in \$(awk 'BEGIN{FS="${params.FIELD_SEP}"} {print \$${params.COL}}' ${ID_file}); do
                        wget -q --header='Content-type:text/x-fasta' "https://rest.ensembl.org/sequence/id/\$i?mask_feature=${params.MASK_INTRON};type=genomic" -O - > \$i.fasta;
                done
                """
        } else if(type_of_input=='T' && type_of_output.contains('cd')) {
                """
                for i in \$(awk 'BEGIN{FS="${params.FIELD_SEP}"} {print \$${params.COL}}' ${ID_file}); do
                        wget -q --header='Content-type:text/x-fasta' "https://rest.ensembl.org/sequence/id/\$i?type=${type_of_output}" -O - > \$i.fasta;
                done
		"""
        } else if(type_of_input=='G' && type_of_output=='protein') {
                """
                for i in \$(awk 'BEGIN{FS="${params.FIELD_SEP}"} {print \$${params.COL}}' ${ID_file}); do
                        wget -q --header='Content-type:text/x-seqxml+xml' "https://rest.ensembl.org/sequence/id/\$i?multiple_sequences=1;type=protein" -O - > \$i.fastaxml;
                done
                """
        } else if(type_of_input=='P') {
	        """
                for i in \$(awk 'BEGIN{FS="${params.FIELD_SEP}"} {print \$${params.COL}}' ${ID_file}); do
                        wget -q --header='Content-type:application/json' "https://rest.ensembl.org/sequence/id/\$i" -O - > \$i.fastaxml;
                done
                """
	}
}



workflow ensembl_types_request_handler {

	take:
	pattern_idfile
	type_in
	type_out

	main:
	
	// error section of missing or wrong inputs
	
	allowed_intypes = ['G', 'P', 'T']
	allowed_outypes = ['genomic', 'cds', 'cdna', 'protein']
	if ( !allowed_intypes.contains(type_in) ) {
		log.info "the --TYPE_INPUT value given is not accepted, value : ${type_in}"
		log.info "allowed : ${allowed_intypes}"
	} else if ( !allowed_outypes.contains(type_out) ) {
		log.info "the --TYPE_OUTPUT value given is not accepted, value : ${type_out}"
                log.info "allowed : ${allowed_outypes}"
	}
	
	idfile = Channel.fromPath(pattern_idfile)
	get_ensembl_emtry(idfile, type_in, type_out)

	emit:
	stout = get_ensembl_emtry.out.standardout
}


workflow {
	ensembl_types_request_handler(params.INPUT_IDS, params.TYPE_INPUT, params.TYPE_OUTPUT)
	ensembl_types_request_handler.out.stout.subscribe  onError: { println it }, onComplete: { println 'Done -> ${params.OUTPUT_DIR}' }
}

