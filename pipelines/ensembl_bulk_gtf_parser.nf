#!/usr/bin/env nextflow



// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline extracts the lines presenting the queried gene ids for a bulk download of ensembl gtf info."
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
	log.info "              on top of that, the requirement for the filenames is also to have the Species_name <- like this at the beginning of the file"
	log.info "              or like Species_name_name and that <optional> does not have a number at the end, this is used to impose the order explained before"
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

params.CONTAINER =  "ubuntu@sha256:2d7ecc9c5e08953d586a6e50c29b91479a48f69ac1ba1f9dc0420d18a728dfc5" // Ubuntu 22.04.1 mawk 1.3.4
params.IN = false
params.FS = '\t'
params.SPECIE_COL = 1
params.GENEID_COL = 2
params.GTF = false
params.OUT_DIR = "${launchDir}/GTF/"
params.OUT_NAME = 'gene'


// Include section
include { awk_pairs } from "${params.PIPES}awk_field_extractor" addParams(CONTAINER: "ubuntu@sha256:2d7ecc9c5e08953d586a6e50c29b91479a48f69ac1ba1f9dc0420d18a728dfc5" ) // the same used in this module




process gtf_lines_extracter  {
	container params.CONTAINER
	tag { "${specie}_${geneID}" }

	input:
	tuple val(specie), val(geneID), path('*')
	val flag_out

	output:
	path "*.gtf", emit: gtf_files
	stdout emit: standardout	//for debug

	script:
	outname = "${geneID}.gtf"
	if ( flag_out=='specie' ) {
		outname = "${specie}.gtf"
	}
	"""
	## forcing to start with species_name.assembly.version.gtf.gz file
	gzip -cd \$(ls *[0-9].gtf.gz) | grep ${geneID} > ${outname} || [[ \$? == 1 ]]		## preventing grep from sending error on not found match 
	if [ -s "${outname}" ]
	then
		exit 0											## exiting with no error
	else
		## listing other remaining files abinitio optional field should be the first lexographically
		for i in \$(ls *.*.*.*.gtf.gz); do gzip -cd \$i | grep 'bubba' > "${geneID}.gtf" || [[ \$? == 1 ]]; done
	fi
	if ! [ -s "${outname}" ]
	then
		rm "${outname}"
		echo "${geneID} has not been found"
	fi
	"""
}


workflow ensembl_gtf_parser  {

	take:
	pattern_to_in
	pattern_to_gtf

	main:

	// error section of missing or wrong inputs

	if ( !pattern_to_in ) {
		log.info "ERROR: no valid input given, pass --IN argument from command line or type --help for description of pipeline"
		exit 1
	} else if ( !pattern_to_gtf ) {
		log.info "ERROR: no valid gtf given, pass --GTF argument from command line or type --help for description of pipeline"
		exit 1
	}

	// Actual pipeline section

	in_files = Channel.fromPath(pattern_to_in)
	awk_pairs(in_files, params.FS, ';', params.SPECIE_COL, params.GENEID_COL, false)
	awk_pairs.out.stout.map{ it -> [(it.split(';')[0]).trim().replace(" ", "_").replace("-", "_"), (it.split(';')[1]).trim()] }.set{ tupled_specie_gene }
	Channel.fromPath(pattern_to_gtf).map{ it -> ["${it.getSimpleName()}", it]}.groupTuple().set{ gtf_files_tuple }
	tupled_specie_gene.join(gtf_files_tuple).set{ correct_matches }
	gtf_lines_extracter(correct_matches, params.OUT_NAME)

	emit:
	stout = gtf_lines_extracter.out.standardout 
	final_out = gtf_lines_extracter.out.gtf_files
}

workflow {
	ensembl_gtf_parser(params.IN, params.GTF)
	ensembl_gtf_parser.out.stout.view()
	ensembl_gtf_parser.out.final_out.view()
}
