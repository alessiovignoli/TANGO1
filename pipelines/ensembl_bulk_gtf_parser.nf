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
        log.info "--OUT_NAME	optional field, default gene, accepted values: <gene> <specie> <specie_nosep> "
	log.info "		tells the script if the name of the output files has to ge the geneID.gtf or Species_name.gtf or Speciesname.gtf (no _)"
        log.info "--CHR		optional field, default false, the extracted lines from the gtf files are going to be only the ones associated with the"
	log.info "		queried geneID. If a value is passed to this flag the extracted line will be instead the whole genes found on the chromosome"
        log.info "		where the queried geneID lies. By Chromosome is also intented scaffold, basically whatever identifier fills the first field"
        log.info "		of the gtf lines of the found queried geneID."
        log.info "		Be aware if this falg is not false the output files will be compressed (gzip), since they tend to be heavy."
	log.info "		And they're anme will change to accomodate the chromosome info, adding it to what is explained above example:"
        log.info "		geneID_chr1.gtf.gz or Species_name_chr1.gtf.gz or Speciesname_chr1.gtf.gz"
        log.info "		The pipeline also outputs a tab separated file with the correspondances species name chromosome of rewuested gen, "
	log.info "		in the same out dir, called specie_chromosome.mapper"
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
params.CHR = false


// Include section
include { awk_pairs } from "${params.PIPES}awk_field_extractor" addParams(CONTAINER: "ubuntu@sha256:2d7ecc9c5e08953d586a6e50c29b91479a48f69ac1ba1f9dc0420d18a728dfc5" ) // the same used in this module




process gtf_lines_extracter  {
	publishDir(params.OUT_DIR, mode: 'move', overwrite: false, saveAs: { filename -> if ( filename.startsWith("TMP")) null
										else filename 
										})
	container params.CONTAINER
	tag { "${specie}_${geneID}" }

	input:
	tuple val(specie), val(geneID), path('*')
	val flag_out
	val chr_switch

	output:
	path "*.gtf*", emit: gtf_files, followLinks: false
	stdout emit: standardout	// for error message and mapper specie chr iD info

	script:
	outname = "${geneID}"
	if ( flag_out=='specie' ) {
		outname = "${specie}"
	} else if ( flag_out=='specie_nosep' ) {
		outname = "${specie}".replace('_', '')
	} else {
		println("the --OUT_NAME value has not been recognized, given: $flag_out{}   allowed : gene, specie, specie_nosep   the output name will be the geneID")
	}

	if ( !chr_switch  ) {
		"""
		## forcing to start with species_name.assembly.version.gtf.gz file
		gzip -cd \$(ls *[0-9].gtf.gz) | grep ${geneID} > ${outname}.gtf || [[ \$? == 1 ]]		## preventing grep from sending error on not found match 
		if [ -s "${outname}.gtf" ]
		then
			exit 0											## exiting with no error
		else
			## listing other remaining files abinitio optional field should be the first lexographically
			for i in \$(ls *.*.*.*.gtf.gz)
			do 
				gzip -cd \$i | grep ${geneID} > ${outname}.gtf || [[ \$? == 1 ]]
				if [ -s "${outname}.gtf" ]
				then
					exit 0
				fi
			done
		fi
		if ! [ -s "${outname}.gtf" ]
		then
			mv ${outname} TMP.gtf
			echo "GENE NOT PRESENT: ${geneID} , in ${specie}.* files"
		fi
		"""
	} else {
		"""
		## First finding which is the first field associated to the geneID
		## forcing to start with species_name.assembly.version.gtf.gz file
		CHR_ID=\$(gzip -cd \$(ls *[0-9].gtf.gz) | grep ${geneID} | head -n 1 | cut -f 1)
		if [ -z \$CHR_ID]								## checking if variable is empty
		then
			## the geneID has not been found in species_name.assembly.version.gtf.gz so the other files will be used for it
			for i in \$(ls *.*.*.*.gtf.gz)
			do
				CHR_ID=\$(gzip -cd \$i | grep ${geneID} | head -n 1 | cut -f 1)
				if ! [ -z \$CHR_ID]
				then
					gzip -cd \$i | grep -P "^\$CHR_ID""\t" > ${outname}.gtf
					gzip  ${outname}.gtf
					echo ${specie} \${CHR_ID}
					exit 0
				fi
			done
		else
			## the geneID has been found in species_name.assembly.version.gtf.gz already 
			gzip -cd \$(ls *[0-9].gtf.gz) | grep -P "^\$CHR_ID""\t" > ${outname}.gtf
			gzip  ${outname}.gtf 
			echo ${specie} \${CHR_ID}
		fi
		if [ -z \$CHR_ID]								## for error message when nothing is found
		then
			touch TMP.gtf								## to avoid error message on missing output
			echo "GENE NOT PRESENT: ${geneID} , in ${specie}.* files"
		fi
		"""
	}
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
	tupled_specie_gene.join(gtf_files_tuple, remainder: true).set{ all_matches } 
	all_matches.filter{ it[1]!=null && it[2]!=null }.set{ correct_matches }
	all_matches.filter{ it[2]==null }.set{ not_found_species }
	//correct_matches.view()
	//not_found_species.view()
	gtf_lines_extracter(correct_matches, params.OUT_NAME, params.CHR)
	gtf_lines_extracter.out.standardout.filter { it.toString().startsWith("GENE NOT PRESENT:") }.set{ stout }
	gtf_lines_extracter.out.standardout.filter { !it.toString().startsWith("GENE NOT PRESENT:") }.collectFile(name: 'specie_chromosome.mapper', storeDir: params.OUT_DIR)

	emit:
	stout 
	final_out  = gtf_lines_extracter.out.gtf_files
	not_found_species
}

workflow {
	ensembl_gtf_parser(params.IN, params.GTF)
	ensembl_gtf_parser.out.stout.view()
	ensembl_gtf_parser.out.not_found_species.subscribe onNext: { println "NOT FOUND SPECIES: ${it[0]} , no file was in the format  ${it[0]}.*" }
}
