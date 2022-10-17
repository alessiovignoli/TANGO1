#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline extracts the lines presenting the queried gene ids for a bulk download of ensembl gtf info."
        log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info '--IN		a file with consistent field separator, where on each line there are at least two column/field'
        log.info "		this two field contain one the species name in scentific format the other an ensemble fasta id code present in such specie"
	log.info "		example:"
        log.info "		first field;	Homo sapiens;	third_field;	1 ;	ecc"
	log.info "              it does not need to be the whole fasta header line, it can be the really unique part of it, (subsequence), aka the fasta identifier"
	log.info "              in the above case is chromosome 1, id -> <1 >, so that the program does not take also 11, 12, 13 ecc... ,it matches with a in python3 statement."
        log.info "--FS		optional field, default <tab>, tells the script how the fields in --IN are separated, the extraction is done with awk"
        log.info "		through the use of FS='' function, look at awk documentation for more details on specifiable patterns."
        log.info "--SPECIE_COL	optional field, default 0 (first column), tells the script which is the column to be extracted as species name."
	log.info "--FASTAID_COL	optional field, default 1 (second column), tells the script which is the column to be extracted as fasta id, as the above flag"
	log.info "		when the number of fields/column is variable look at the usage of NF value in awk for referring to last field. pass NF as value for last column."
	log.info "--FASTA 	path to the file/s in fasta format, they have to be zipped in .gz extension otherwise nothing will be extracted."
        log.info "		the fasta files downloadable from ensembl follow this naming:"
        log.info "		<species>.<assembly>.<sequence type>.<id type>.<id>.gtf.gz"
        log.info "		<species>:		The systematic name of the species."
        log.info "		<assembly>:		The assembly build name."
        log.info "		<sequence type>:	The type of masking done on the DNA."
	log.info "		<id type>:		optional, variable field representing the level of the assembly, chromosome seqlevel ecc.."
	log.info "              <id>:			optional, The actual sequence identifier. Depending on the <id type>."
	log.info "              for more details look at ->  http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/README"
        log.info "		the files given to this pipeline with this flag are expected to folow the same rules."
	log.info "              on top of that, the requirement for the filenames is also to have the Species_name <- like this in <species>"
	log.info "              (should be default in Ensembl) or like Species_name_name for subspecies."
	log.info "              If more than one file is found to match this rules for the same species name, the searched file will be the one with no masking:"
	log.info "              <sequence type> = dna  , this is also searched by the pipeline and if not found an error message missing file will prompt."
        log.info "--OUT_DIR	optional field, default launchDir/Genomes, launchDir is a nextflow variable representing where the pipeline has been launched from."
        log.info "--OUT_NAME	optional field, default specie, accepted values: <specie> <specie_nosep>  "
	log.info "		tells the script if the name of the output files has to be the Species_name_<fastaID>.fa.gz or Speciesname_<fastaID>.fa.gz (no _)"
        log.info ""
        log.info ""
        log.info ""
        log.info ""
        exit 1
}

params.CONTAINER =  "ubuntu@sha256:2d7ecc9c5e08953d586a6e50c29b91479a48f69ac1ba1f9dc0420d18a728dfc5" // Ubuntu 22.04.1 mawk 1.3.4
params.IN = false
params.FS = '\t'
params.SPECIE_COL = 1
params.FASTAID_COL = 2
params.FASTA = false
params.OUT_DIR = "${launchDir}/Genomes/"
params.OUT_NAME = 'specie'


// Include section
include { awk_pairs } from "${params.PIPES}awk_field_extractor" addParams(CONTAINER: "ubuntu@sha256:2d7ecc9c5e08953d586a6e50c29b91479a48f69ac1ba1f9dc0420d18a728dfc5" ) // the same used in this module




workflow ensembl_genome_parser  {

	take:
	pattern_to_in
	pattern_to_fa

	main:

	// error section of missing or wrong inputs

	if ( !pattern_to_in ) {
		log.info "ERROR: no valid input given, pass --IN argument from command line or type --help for description of pipeline"
		exit 1
	} else if ( !pattern_to_fa ) {
		log.info "ERROR: no valid gtf given, pass --FASTA argument from command line or type --help for description of pipeline"
		exit 1
	}

	// Actual pipeline section

	in_files = Channel.fromPath(pattern_to_in)
	awk_pairs(in_files, params.FS, ';', params.SPECIE_COL, params.FASTAID_COL, false)
	awk_pairs.out.stout.map{ it -> [(it.split(';')[0]).trim().replace(" ", "_").replace("-", "_"), (it.split(';')[1]).replace("\n", "")] }.set{ tupled_specie_fastaID }
	Channel.fromPath(pattern_to_fa).map{ it -> ["${it.getSimpleName()}", it]}.groupTuple().set{ fasta_files_tuple }
	tupled_specie_fastaID.join(fasta_files_tuple, remainder: true).set{ all_matches } 
	all_matches.filter{ it[1]!=null && it[2]!=null }.set{ correct_matches }
	all_matches.filter{ it[2]==null }.set{ not_found_species }
	correct_matches.view()
	//not_found_species.view()
	//gtf_lines_extracter(correct_matches, params.OUT_NAME, params.CHR)

	emit:
	stout = 'bubba' //gtf_lines_extracter.out.standardout 
	//final_out = gtf_lines_extracter.out.gtf_files
	not_found_species
}

workflow {
	ensembl_genome_parser(params.IN, params.FASTA)
	ensembl_genome_parser.out.stout.view()
	ensembl_genome_parser.out.not_found_species.subscribe onNext: { println "NOT FOUND SPECIES: ${it[0]} , no file was in the format  ${it[0]}.*" }
}
