#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "this pipeline computes the number of entities that belong to each kingfom of life"
        log.info ""
        log.info "the input --INPUT is either a fasta file or a file with the following specfics:"
        log.info "	>prot_id taxid=9606"
        log.info "so the hwader has to start with a > sign as first charachter and in the same line there has to be "
        log.info "a taxid field, like uniptot headers have, the file is going to look for a list of different taxid field "
        log.info "keywords, see        ncbi_file_querier.py     for more details"
        log.info "once the script identifies the taxid keyword it splits the line on it and retrieves the actual number digit by digit"
        log.info "this should allow some flexibility with the input "
	log.info "the other mandatory field is the database the type accepted is the one of NCBI categories.dmp taxonomy where the strusture is what follows: (first tab is just for visual clarity)"
        log.info "	B	251701  323"
        log.info "	B	329     329"
        log.info "the first field is the letter key for the kingdom, the second the sub-specie taxid and the third the species taxid"
        log.info "--NCBI_DB flag has to b e used to give this file"
        log.info "if some taxid are not given then the process either otputs the results and the unfound taxid"
        log.info "or it searches the unfound on the other database given see below"
        log.info ""
        log.info "the other optional field is the database the type accepted is the one of NCBI taxonomy nodes.dmp where the strusture is what follows: (first tab is just for visual clarity)"
        log.info "	147553	|	451866	|	class	|		|	4	|	1	|	1	|	1	|	4	|	1|	0	|	0	|		|"
        log.info ""
        log.info ""
        log.info "--NCBI_FULL_DB flag has to b e used to give this file"
        log.info "it can be downloaded from https://ftp.ncbi.nih.gov/pub/taxonomy/"
        log.info "and it structure is explained here:"
        log.info "tax_id		--1 field	-- node id in GenBank taxonomy database"
        log.info "parent tax_id		--2 field	-- parent node id in GenBank taxonomy database"
        log.info "rank			--3 field	-- rank of this node (superkingdom, kingdom, ...)"
        log.info "embl code		--4 field	-- locus-name prefix; not unique"
        log.info "division id		--5 field	-- the field the script is interested in"
        log.info "for the other fieldslook at the readme file in the above link"
        log.info ""
        log.info "the values of the division id field are the following and interpreted as described:"
	log.info "0	|	BCT	|	Bacteria	|	Bacteria	|"
        log.info "1	|	INV	|	Invertebrates	|	Eukaryota	|"
        log.info "2	|	MAM	|	Mammals	|	Eukaryota	|"
        log.info "3	|	PHG	|	Phages	|	Viruses	|"
        log.info "4	|	PLN	|	Plants and Fungi	|	Eukaryota 	|"
        log.info "5	|	PRI	|	Primates	|	Eukaryota	|"
        log.info "6	|	ROD	|	Rodents	|	Eukaryota	|"
        log.info "7	|	SYN	|	Synthetic and Chimeric	|	Eukaryota 	|"
        log.info "8	|	UNA	|	Unassigned	|	Unassigned	|"
        log.info "9	|	VRL	|	Viruses	|	 Viruses	|"
        log.info "10	|	VRT	|	Vertebrates	|	Eukaryota 	|"
        log.info "11	|	ENV	|	Environmental samples	|	Anonymous sequences cloned directly from the environment	|"
        log.info ""
        log.info "### 	WARNING  ####    archea are classified as Bacterias in nodes.dmp"
        exit 1
}

params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.NCBI_DB = "${params.TEST_DIR}bubbabubba" 
params.NCBI_FULL_DB = false

process ncbi_searcher {
    //publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
    container params.CONTAINER
    tag { "${in_fasta}" }

    input:
    path in_fasta 
    path in_ncbi
    path pyscript

    output:
    stdout emit: standardout            // used only for check during debugging

    script:
    """
    ./${pyscript} ${in_fasta} ${in_ncbi}
    """
}

workflow ncbi_taxa_parser {

    take:
    pattern_to_input
    pattern_to_db

    main:
    in_id = Channel.fromPath(pattern_to_input)
    in_db = Channel.fromPath(pattern_to_db)
    fastalike_ncbi_search = params.SCRIPTS + "ncbi_file_querier.py"
    ncbi_searcher(in_id, in_db, fastalike_ncbi_search)

    emit:
    stout = ncbi_searcher.out.standardout
}


workflow {
    if ( params.NCBI_DB != false ) {
        ncbi_taxa_parser(params.INPUT, params.NCBI_DB)
        ncbi_taxa_parser.out.stout.view()
    }
}
