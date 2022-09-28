#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline ecxtract the lines presenting the queried gene ids for a bulk download of ensembl gtf info."
        log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info '--IN		a file with consistent field separator, where on each line there are at least two column/field'
        log.info "		this two field contain one the species name in scentific format the other an ensemble gene id code present in such specie"
	log.info "		example:"
        log.info "		first field;	Homo sapiens;	third_field;	ENSG00000279493;	ecc"
        log.info "--FS		optional field, default <tab>, tells the script how the fields in --IN are separated, the extraction is done with awk"
        log.info "		through the use of FS='' function, look at awk documentation for more details on specifiable patterns."
        log.info "--SPECIE_COL	optional field, default 0 (first column), tells the script which is the column to be extracted as species name."
	log.info "--GENEID_COL	optional field, default 1 (second column), tells the script which is the column to be extracted as gene id, as the above flag"
	log.info "		when the number of fields/column is variable look at the usage of NF value in awk for referring to last field. pass NF as value for last column."
	log.info "--GTF		path to the file/s in gtf format, they have to be zipped in .gz extension otherwise nothing will be extracted."
        log.info "		the gtf files downloadable from ensembl follow this naming:"
        log.info "		<species>.<assembly>.<version>.<optional>.gtf.gz"
        log.info "		<species>:       The systematic name of the species."
        log.info "		<assembly>:      The assembly build name."
        log.info "		<version>:       The version of Ensembl from which the data was exported."
	log.info "		<optional>:      The additional flag, usually abinitio, can be other read ensembl documentation for more explanation."
        log.info "		the files given to this pipeline with this flag are expected to folow the same rules."
	log.info "		when more than one file is found to start with species name the hierachy is the following:"
	log.info "		without the optional field -> with abinitio as optional field -> all other optional fields in sorting order."
	log.info "              Once a geneID is found in a file the search is stopped for that gene, aka the other files in lower hirerchy are not searched."
        log.info "--OUT_DIR	optional field, default launchDir/GTF, launchDir is a nextflow variable representing where the pipeline has been launched from."
        log.info "--OUT_NAME	optional field, default gene, accepted values: <gene> <specie>  "
	log.info "		tells the script if the name of the output files has tio ge the geneID.gtf or Species_name.gtf"
        log.info "####   WARNING #####"
        log.info "the output are as many files as there are queried species, and they're name is the following:"
        log.info "<species>.<assembly>.<version>."
	log.info ""
        log.info ""
        log.info ""
        log.info ""
	log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "alessiovignoli3/tango-project@sha256:57013bf372b519245608c95fd60a38f9e5d65775aaf18c2d711031818c1a145e" // awk and wget
params.IN = false
params.FS = '\t'
params.SPECIE_COL = 0
params.GENEID_COL = 1
params.GTF = false
params.OUT_DIR = "${launchDir}/GTF/"
params.OUT_NAME = 'gene'

/*

process fasta_longest_picker  {
	container params.CONTAINER
	tag { "${fasta}" }

	input:
	path fasta
	path pyscript

	output:
	path "*.tmp", emit: tmp_files
	//stdout emit: standardout	for debug

	script:
	"""
	python3 ${pyscript} --fasta ${fasta} > ${fasta}.tmp
	"""
}
*/

workflow ensembl_gtf_parser  {

	take:
	pattern_to_in
	pattern_to_gtf

	main:

	// error section of missing or wrong inputs

	if ( !pattern_to_in ) {
		log.info "no  --help for description of pipeline"
	}


	//emit:
}

workflow {
	ensembl_gtf_parser(params.IN, params.GTF)
}
