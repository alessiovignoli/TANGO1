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
        log.info "keywords, see        # # # #  .py     for more details"
        log.info "once the script identifies the taxid keyword it splits the line on it and retrieves the actual number digit by digit"
        log.info "this should allow some flexibility with the input "
        log.info ""
        log.info "the other optional field is the database the type accepted is the one of NCBI taxonomy where the strusture is what follows: (first tab is just for visual clarity)"
        log.info "	B	251701  323"
        log.info "	B	329     329"
        log.info "the first field is the letter key for the kingdom, the second the sub-specie taxid and the third the species taxid"
        log.info "--NCBI_DB flag has to b e used to give thid file"
        log.info "if no database file is given or if taxid retrieved from the input file are not found in the dbfile"
        log.info "they are searched against ena rest api that is again based on NCBI database"
        log.info "if there are still ufound taxid they are printed to a file <input filename>_notfound.taxids"
        log.info "in the directory where the script is lounched"
        exit 1
}

params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.NCBI_DB = false

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
