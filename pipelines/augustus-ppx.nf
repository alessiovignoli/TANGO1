#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline execute the uagustus ppx function, basically a prot-MSA refined gene prediction"
        log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info '--IN		file or filepath as glob pattern to genome files containing portion of DNA sequence in fasta format'
        log.info "		when many files are given, each one is treated as separate genome/specie."
	log.info "--HS		optional flag, default <space>, tells the script which is the field separator oh the header in genome file"
	log.info "--ID_POS	optinal flag, default 1, first field, toghther with the above flags specifies which is the ID present in the"
	log.info "		header of the genome file. This is necessary for the retrival of the species name associated to the genome file"
	log.info "		same concepts appliable for --FS and --POS hold. read below."
	log.info "--SPECIES	file containing on each line some fields separated by the same separator, among the fields there has to be"
        log.info "		a scientific species name like    Homo sapiens    and (on the same line) the id present in the header"
        log.info "		of the corresponding genome file, exxample:"
        log.info "		genome file   ->    homo_chr1.fa    content:"
        log.info "		>IDcode1 some miscellanues info"
        log.info "		AACATGGCTGCGGCGCCTGGGCTGCTAGTCTG....."
        log.info "		species file   ->   species.csv	    content:"
        log.info "		IDcode1,something,something,Homo sapiens"
        log.info "		IDcode2,some,thing,Gorilla gorilla"
	log.info "		..."
        log.info ""
        log.info "		if more than one species file is given all will be searched but the program expetcts"
	log.info " 		them to have the same format, most important same field separator and same position of the scientific name."
	log.info "              Only the first match is retained."
        log.info "--FS		optional flag, default <tab>, flag used to tell the programm what is the field separator "
	log.info "		used in the species file, to give something different pass --FS ';\\t' ,"
        log.info "		in this example the columns are eparated by a comma followed by a <tab>, quotation are necessary in this case."
        log.info "		The extraction of the field is one through awk program, specifically FS= variable, look at awk documentation."
        log.info "--POS		optional flag, default 1, first field, it tells the script which is the field position to extract the species name"
	log.info "		last field can be specified as NF and second to last as NF-1, llok at awk documentation"
        log.info "--REF_FA	give this flag if the sequences comprising the MSA used for enhanche the genome predicitons have to be"
	log.info "		alligned still. This flags accepts only one fasta file, so if sequences nedded for alignment are in differnt"
        log.info "		files is uop the user to put them in one. The produced alignment is going to be published in the OUT_DIR directory"
        log.info "--REF_ALN	pass this flag if the sequences are already alligned, the format accepted by augustus are: FASTA or CLUSTAL."
	log.info "              This flag accepts only one file."
	log.info "		"
	log.info "the above two flags are mutually exclusive, only one of them can be given."
        log.info "##########   WARNING    ###########"
	log.info "if more than one file is passed in one of the above --REF flags, for example as a glob path"
	log.info "only the first file is going to be used"
        log.info ""
        log.info "--PREP_ALN	optional flag, default false no neeed for   prepareAlign   module of Augustus."
        log.info "		basically if the input MSA or the one generated is a bit too gappy or messy or big"
	log.info "		then is better to give this flag and use such module specifically made for circumventing the problem."
        log.info "		More info can be found at    https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md"
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "alessiovignoli3/tango-project@sha256:9a351679d2f41c54b2baabb44feb415c870344406f11bde627854d97f810aaf9" // augustus 3.4.0
params.IN = false
params.HS = ' '
params.ID_POS = 1
params.SPECIES = false
params.FS = '\t'
params.POS = 1
params.REF_FA = false
params.REF_ALN = false
params.OUT_DIR = "${params.TEST_DIR}augustus_results/"
params.PREP_ALN = false

// Modules dependencie section

include { oneliner_ch } from "${params.PIPES}fasta_oneliner" addParams(PUBLISH: "false", CONTAINER: "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43")	// the same used in this module
include { align_generation } from "${params.PIPES}aln_generation_colouring" addParams(OUTPUT_DIR: params.OUT_DIR) 
include { clustal_converter } from "${params.PIPES}from_clustal_to_fasta_aln" addParams(PUBLISH: "false", CONTAINER: "cbcrg/tcoffee@sha256:36c526915a898d5c15ede89bbc3854c0a66cef22db86285c244b89cad40fb855")  // what is necessary for such module


process  aug_prep_aln {
	container params.CONTAINER
	tag { "${aln}" }

	input:
	path aln

	output:
	//path "*", emit: tmp_files
	stdout emit: standardout	// for debug

	script:
	outname = "${aln}".split('\\.')[0] + '.preppedaln'
	"""
	prepareAlign < ${aln} > ${outname} 
	"""
}


process aug_msa_to_profile {
	container params.CONTAINER
        tag { "${msa}" }

        input:
        path msa

	output:
        path "*.prfl", emit: profile_file
        stdout emit: standardout        // for debug

        script:
        outname = "${msa}".split('\\.')[0] + '.prfl'
        """
	msa2prfl.pl ${msa} > ${outname}
        """
}



process aug_ppx {
	publishDir(params.OUT_DIR, mode: 'move', overwrite: false)
	container params.CONTAINER
        tag { "${dna}" }

	input:
	path dna
	path profile
	path species_file

	output:
	//path "*", emit: 
        stdout emit: standardout        // for debug

	script:
	outname = "${dna}".split('\\.')[0] + '.'
	"""
	echo ${aug_species_list}
	#ID_HEADER_GENOME=\$(grep '>' ${dna} | head -n 1 | cut -d '>' -f 2| awk 'BEGIN{FS="${params.HS}"} {print \$${params.ID_POS}}') 
	#SPECIES_NAME=\$(grep "\$ID_HEADER_GENOME" ${species_file} | awk 'BEGIN{FS="${params.FS}"} {print \$${params.POS}}')
	#grep "\$SPECIES_NAME" tmp
	#augustus --proteinprofile=${profile} --species=\$SPECIES_NAME ${dna}
	"""
}




workflow  augustus_ppx {

	take:
	pattern_to_genomes
	pattern_to_species
	pattern_to_fasta
	pattern_to_msa

	main:

	// error section of missing or wrong inputs

	if ( !pattern_to_genomes ) {
		log.info "ERROR: --IN is a mandatory flag, no file specified, for more details type --help"
                exit 1
	} else if ( !pattern_to_species ) {
		log.info "ERROR: --SPECIES is a mandatory flag, no file specified, for more details type --help"
		exit 1
	} else if ( !pattern_to_fasta && !pattern_to_msa ) {
		log.info "ERROR: nor --REF_FA nor --REF_ALN has been given, one of the two is mandatory, for more details type --help"
		exit 1
	} 

	
	// alignment preparation section

	in_ref = null
	if ( pattern_to_fasta ) {
		fastas = channel.fromPath(pattern_to_fasta)
		oneliner_ch(fastas.first(), 'onelinefa')
		align_generation(oneliner_ch.out.onelinefasta)
		in_ref = align_generation.out.aln_file
	} else {
		in_ref = channel.fromPath(pattern_to_msa)
	}
	prepped_aln = in_ref
	if ( params.PREP_ALN ) {
		clustal_converter(in_ref, '.alnfasta')
		prepped_aln = aug_prep_aln(clustal_converter.out.fasta_aln)
	}

	// actual augustus prediction section

	aug_msa_to_profile(prepped_aln)
	genome = channel.fromPath(pattern_to_genomes)
	species_id_match = channel.fromPath(pattern_to_species)
	aug_ppx(genome, aug_msa_to_profile.out.profile_file, species_id_match)

	emit:
	stout = aug_ppx.out.standardout
}

workflow {
	augustus_ppx(params.IN, params.SPECIES, params.REF_FA, params.REF_ALN)
	augustus_ppx.out.stout.view()
}

