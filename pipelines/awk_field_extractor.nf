#!/usr/bin/env nextflow

// ###
// ###
// ###  IMPORT named workflows not    awk_mode_selector
// ###
// ###


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline extracts the field requensted from the input file using awk, it is not really thought to be launched by itself. It is"
        log.info "written as a collection of importable workflows, Having said that is also possible to use it by itself from command line."
	log.info ""
        log.info "Here is the list of flags accepted by the pipeline:"
        log.info '--IN		a file with consistent field separator, where on each line there are at least the number of requested column/field'
        log.info "--FS		optional field, default <tab>, tells the script how the fields in --IN are separated, the extraction is done with awk"
        log.info "		through the use of FS='' function, look at awk documentation for more details on specifiable patterns."
	log.info "--OFS		optional field, default <tab>, tells the script how the otput fields printed by awk have to be separated"
	log.info "              through the use of OFS='' function, look at awk documentation for more details on specifiable patterns."
	log.info "--MODE	optional field, default <single>, tell the script how to use awk, accepted values:"
        log.info "		single, couple, triple, slice"
        log.info "		single:		extracts only the desired column --COL1 from --IN file for every line				## NOT IMPLEMENTED YET"
	log.info "		couple:		extracts only the two column requested --COL1, --COL2 from --IN file for every line"
        log.info "		triple:		extracts only the 3 column requested --COL1, --COL2, --COL3 from --IN file for every line	## NOT IMPLEMENTED YET"
        log.info "		slice:		extracts the n columns specified by the       flag  from --IN file for every line		## NOT IMPLEMENTED YET"
        log.info ""
        log.info "--COL1	optional field, default 1 (first column), tells the script which is the column to be extracted as first argument of output."
	log.info "--COL2	optional field, default 2 (second column), tells the script which is the column to be extracted as second argument of output."
        log.info "--COL3	optional field, default 3 (third column), tells the script which is the column to be extracted as third argument of output."
        log.info "		when the number of fields/column is variable look at the usage of NF value in awk for referring to last field. pass NF as value for last column."
	log.info "--OUT_TYPE	optional field, default <false>, tells the script if the output has to be a standard out nextflow channel, or a file channel."
        log.info "		in case an argument is passed that is going to be the suffix of the output file as: <--IN>_<--OUT_TYPE>.txt "
        log.info ""
        log.info ""
	log.info ""
        log.info ""
        exit 1
}

params.CONTAINER = "ubuntu@sha256:2d7ecc9c5e08953d586a6e50c29b91479a48f69ac1ba1f9dc0420d18a728dfc5" // Ubuntu 22.04.1 mawk 1.3.4
params.IN = false
params.FS = '\t'
params.OFS = '\t'
params.MODE = 'single'
params.COL1 = 1
params.COL2 = 2
params.COL3 = 3
params.OUT_TYPE = false



process awk_pair_extracter  {
	container params.CONTAINER
	tag { "${infile}" }

	input:
	path infile
	val field_separator
	val out_field_separator
	val col1_pos
	val col2_pos
	val suffix
	
	output:
	path "*.txt", emit: awk_file
	stdout emit: standardout	// for debug

	script:
	outname = "${infile.getSimpleName()}" + '_' + "${suffix}.txt"
	"""
	awk 'BEGIN{FS="${field_separator}"; OFS="${out_field_separator}"} {print \$${col1_pos}, \$${col2_pos}}' ${infile} > ${outname}
	"""
}


// ###
// ###
// ###			STUFF TO IMPORT HERE BELOW  ||
// ###			                            \/
// ###


workflow awk_pairs {
	
	take:
	channeled_in_files
	field_sep
	out_field_sep
	first_col_pos
	second_col_pos
	out_flag

	main:
	awk_pair_extracter(channeled_in_files, field_sep, out_field_sep, first_col_pos, second_col_pos, out_flag)
	
	// handling of output
	stout = awk_pair_extracter.out.awk_file
	if ( out_flag==false || out_flag=='false' || out_flag=='False') {
		stout = awk_pair_extracter.out.awk_file.splitText()
	}

	emit:
	stout 
}



// ###
// ###			                            /\
// ###                  STUFF TO IMPORT HERE ABOVE  ||
// ###                                             
// ###





workflow awk_mode_selector {

        take:
	pattern_to_in
	mode_specifier
	out_type_specifier

	main:
	
	// error section of missing or wrong inputs

	allowed_modes = ['single', 'couple', 'triple', 'slice']		// check help section to see if implemented yet
	if ( !pattern_to_in ) {
		log.info "ERROR: no valid input given, pass --IN argument from command line or type --help for description of pipeline"
		exit 1
	} else if ( !allowed_modes.contains(mode_specifier) ) {
		log.info "ERROR: given mode specifier not recognized, --MODE value: ${mode_specifier}"
		log.info "ERROR: allowed values: ${allowed_modes}"
		log.info 'for more details type --help or take a look at conf/pd.config'
                exit 1
	}

	// Part of the various modules that have to be imported but called here to allow CL run of pipeline
	
	in_files = Channel.fromPath(pattern_to_in)
	out_awked = ''
	if ( mode_specifier=='single' ) {
		out_awked  = 'single mode not yet implemented   sorry'		// to change
	} else if (  mode_specifier=='couple' ) {
		awk_pairs(in_files, params.FS, params.OFS, params.COL1, params.COL2, out_type_specifier)
		out_awked = awk_pairs.out.stout
	} else if (  mode_specifier=='triple' ) {
		out_awked  = 'triple mode not yet implemented   sorry'          // to change
	} else {
		out_awked  = 'slice mode not yet implemented   sorry'          // to change
	}

	emit:
	out_awked
}

workflow {
	awk_mode_selector(params.IN, params.MODE, params.OUT_TYPE)
	awk_mode_selector.out.out_awked.view()
}
