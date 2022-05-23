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
        log.info "7	|	SYN	|	Synthetic and Chimeric	|	other 	|"
        log.info "8	|	UNA	|	Unassigned	|	Unassigned	|"
        log.info "9	|	VRL	|	Viruses	|	 Viruses	|"
        log.info "10	|	VRT	|	Vertebrates	|	Eukaryota 	|"
        log.info "11	|	ENV	|	Environmental samples	|		|"
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
    stdout emit: standardout           
    path "${outname_404}", emit: not_found_taxids

    script:
    outname_404 = "not_found_in_species-" + "${in_fasta}".split('\\.')[0] + ".headers" 		// DO NOT CHANGE THIS USED LATER AS ID
    """
    ./${pyscript} ${in_fasta} ${in_ncbi} ${outname_404}
    """
}

process ncbi_nodes_searcher {
    //publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
    container params.CONTAINER
    tag { "${in_fastalike}" }

    input:
    path in_fastalike
    path in_ncbi_db
    //val already_computed_kingdoms
    path pyscript

    output:
    stdout emit: standardout            
    path "${out_404}", emit: not_found_taxids2

    script:
    tmp = "${in_fastalike}".split('-').size()
    out_404 = "not_found_in_nodes-"					// DO NOT CHANGE THIS
    for(int i = 1;i<tmp;i++) {						// IT iIS uSED LATER AS IDENTIFIER
        out_404 += "${in_fastalike}".split('-')[i]
    }
    """
    ./${pyscript} ${in_fastalike} ${in_ncbi_db} ${out_404} "true"
    """
}


process ena_rest_api_xml {
    //publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
    container 'alessiovignoli3/tango-project:bash_curl@sha256:82126097c7ed5f374a7e0c0ca4dcd92f7a7678bfac35aa60aa52abcfd9443401'
    tag { "${not_found_ids}" }

    input:
    path not_found_ids

    output:
    stdout emit: standardout

    script:
    """
    Bacteria=0
    Archaea=0
    Eukaryota=0
    Virus=0
    Unclassified=0
    for i in `cat "${not_found_ids}"`; do STR=`curl https://www.ebi.ac.uk/ena/browser/api/xml/\$i | grep 'rank="superkingdom"'`; 
    if [[ \$STR == *"2759"* ]]; then
        Eukaryota=\$((\$Eukaryota + 1))
    elif [[ \$STR == *"10239"* ]]; then
        Virus=\$((\$Virus + 1))
    elif [[ \$STR == *"2157"* ]]; then
        Archaea=\$((\$Archaea + 1))
    elif [[ \$STR == *"2"* ]]; then
        Bacteria=\$((\$Bacteria + 1))
    else
       Unclassified=\$((\$Unclassified + 1))
    fi; done
    echo  "${not_found_ids}"
    echo ["'B'", \$Bacteria, "'A'", \$Archaea, "'E'", \$Eukaryota, "'V'", \$Virus, "'U'", \$Unclassified, "'O'", 0]
    echo 'total taxid found =   '`expr \$Bacteria + \$Archaea + \$Eukaryota + \$Virus`
    """

}

process stdout_collecter {
    //publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
    container params.CONTAINER
    tag { "${suffix_outname}" }

    input:
    val all_stdout
    val suffix_outname

    output:
    path "${suffix_outname}_kingdom_percentages.txt", emit: final_out
    //stdout emit: final_out

    script:
    """
    #!/usr/bin/env python3 

    tmp_list = "${all_stdout}".split('], [')
    #print(tmp_list[0:2])
    l1 = ['B', 0, 'A', 0, 'E', 0, 'V', 0, 'U', 0, 'O', 0]
    l2 = ['B', 0, 'A', 0, 'E', 0, 'V', 0, 'U', 0, 'O', 0]
    l3 = ['B', 0, 'A', 0, 'E', 0, 'V', 0, 'U', 0, 'O', 0]
    overall = ['B', 0, 'A', 0, 'E', 0, 'V', 0, 'U', 0, 'O', 0]
    tot = 0
    for elem in tmp_list:
        actual_list = (elem.split(', [')[1].split('], ')[0]).split(', ')
        id = elem.split(', [')[0]
        tot += int(elem.split('=  ')[1].split(']')[0])
        if 'not_found_in_species' in id:
            for n in range(1, len(actual_list), 2):
                l2[n] += int(actual_list[n])
                overall[n] += int(actual_list[n])
        elif 'not_found_in_nodes' in id:
            for n in range(1, len(actual_list), 2):
                l3[n] += int(actual_list[n])
                overall[n] += int(actual_list[n])
        else:
            for n in range(1, len(actual_list), 2):
                l1[n] += int(actual_list[n])
                overall[n] += int(actual_list[n])
    with open("${suffix_outname}_kingdom_percentages.txt", 'w') as outfile:
        outfile.write('ncbi species search results :\\n' + str(l1) + '\\nncbi nodes search results :\\n' + str(l2) + '\\neni rest api search results :\\n' + str(l3) + '\\noverall results :\\n' + str(overall) + '\\ntotal sequences found :' + str(tot))
    """
}

workflow ncbi_taxa_parser {

    take:
    pattern_to_input
    pattern_to_db
    pattern_nodes_db

    main:
    in_id = Channel.fromPath(pattern_to_input)
    in_db = Channel.fromPath(pattern_to_db)
    fastalike_ncbi_search = params.SCRIPTS + "ncbi_file_querier.py"
    ncbi_searcher(in_id, in_db.first(), fastalike_ncbi_search)
    tmp = ncbi_searcher.out.standardout.collect{ it.split('\\n') }
    stout = false 								
    not_found = false
    final_out = false
    if ( pattern_nodes_db != false ) {
        in_nodes_db = Channel.fromPath(pattern_nodes_db)
        ncbi_nodes_searcher(ncbi_searcher.out.not_found_taxids, in_nodes_db.first(), fastalike_ncbi_search)
        tmp1 = ncbi_nodes_searcher.out.standardout.collect{ it.split('\\n') }
        //tmp1.view()
        not_found = ncbi_nodes_searcher.out.not_found_taxids2
        //ena_rest_api_xml(ncbi_nodes_searcher.out.not_found_taxids2)
        //tmp2 = ena_rest_api_xml.out.standardout.collect{ it.split('\\n') }
        //tmp2.view()
        stout = tmp.concat( tmp1 ).collect()
        //stout.view()
    } else {
        stout = tmp
        not_found = ncbi_searcher.out.not_found_taxids
    }
    prefix = pattern_to_input.split('/').size()
    suffix = 'bubba'
    if ( prefix == 1) {
        suffix = pattern_to_input.split('\\*')[0]
    } else {
        suffix = pattern_to_input.split('/')[-1].split('\\*')[0]
    }
    stdout_collecter(stout, suffix)
    final_out = stdout_collecter.out.final_out

    emit:
    //stout
    not_found
    final_out
}

workflow {
    ncbi_taxa_parser(params.INPUT, params.NCBI_DB, params.NCBI_FULL_DB)
    //ncbi_taxa_parser.out.stout.view()
    ncbi_taxa_parser.out.not_found.view()
    ncbi_taxa_parser.out.final_out.view()
}

