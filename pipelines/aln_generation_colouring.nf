#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'Here is a brief description of the pipeline itself  : '
        log.info '	Here there is the collection of the flags _ parameters specifiable from command line'
        log.info '	if the --INPUT is a collection of fasta headers use the flag --SEARCHFILES followed by filenames or paths separated by a comma'
        log.info '	these are the files on which the headers are compared to retrieve the associated sequences for example:'
        log.info '	    --SEARCHFILES "/home/alessio/Desktop/erasmus-internship/backups/all_mia3_got_so_far_oct14.fasta,/home/alessio/Desktop/erasmus-internship/backups/only_Vertebrata.fasta"'
        log.info '	only complete path accepted, for more info look at the script <from_header_to_fasta.py>'
        log.info '	otherwise --INPUT has to be a fasta file with the sequences needed for the allignment'
        log.info ' 	REMEMBER this pipeline works with one input file at the time for two inputs launch it two times'
        log.info ' '
        log.info '	if the above flag is not given then the pipeline assumes that the input is a fasta file with sequences on one line' 
        log.info '	The signal peptide parameter is mandatory if --TRIMM option is specified'
        log.info '      --TRIMM this can be an integer or two intger separated by a comma whre the first one is cut on the left and the second on the right '
        log.info '      like --TRIMM 50,20'
        log.info '      --SIGNALPEPT this is an integers it has to go with trimm parameter, it tells the script to treat as signalpeptide '
        log.info '      the first n residue of every sequence, it is to avoid to trimm the sequence right at the beginning on a signalpept'
        log.info '      recognized by phobius as TM'
        log.info '      --TRIMM_PLP_DIR is a param utilized in the trimm step to know where to go and look at plp files to use as reference'
        log.info '      for the cut, it can be specified by command line or check the nextflow.config file '
        log.info ''
        log.info '	######### WARNING ############'
        log.info '      --PLP_DIR is a mandatory argument for this script and is the directory where all the plp files of the'
        log.info '	 input sequences are, all the sequences we want to allign have to have an pre-computed plp file associated to them'
        log.info '      it is The folder containning all the results obtained from phobious for the sequences that we want to allign and colour'
        log.info '      The sequences are coloured on the basis of the postirior probability file generated from phobius (.plp)'
        log.info '      This plp files have to be computed before this script'
        log.info '      '
        log.info '      --SPECIAL_HELIX optional parameter used in colouring phase for differential colouring, it can be true or false'
        log.info '	default false, no coloring for special helix'
        log.info '	--MASK is an optional field that tells the script from which position on the sequence has to start the colouring'
        log.info '	 and where to finish, specified extemity are coloured, this is a range made of two integers separated from a comma'
        log.info '	--HYDRO_SCALE to specify which hydrophobicity scale to use in the colouring step'
        log.info '	the supported arguments are kyte for Kyte-Doolittle, GES for GES scale, '
        log.info '	UHS for Unified Hydrophobicity Scale (default),  and <majority> for the majority scale that i invented ;) '
        log.info '	--PHOB_STDOUT this flag is mandatory since the script has to know where to look for the colouring of'
        log.info '	alignment based on phobius stdoutput default format (-long).'
        exit 1
}


params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // slim buster
params.INPUT = "${params.TEST_DIR}bubbabubba"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.TRIMM = false
params.SIGNALPEPT = 0
params.MASK = false
params.SPECIAL_HELIX = false
params.CREATE_SH_REF = false
params.SEARCHFILES = false
params.PLP_DIR = false
params.PHOB_STDOUT = false
params.HYDRO_SCALE = "UHS"

// Modules dependencie section

include { oneliner_ch } from "${params.PIPES}fasta_oneliner" addParams(PUBLISH: "false")
include { converger } from "${params.PIPES}output_files_uniter" 

/*
#phobiuos_results_folder = "${params.parent_ofall_dir}all_phobius_results/MIA2_vertebrate_orthologs_OM/"    // to change if needed remeber to put the slash at the end
#phobius_stdout_filepath = phobiuos_results_folder + ""   // to change , this file is created from the simple re-direction of phobious stdout default, is neede in the phob_stout_colours section. The set originating this prediction has to be the set of sequences presnt in <input_header_file> input file
#params.mask = false            // this is a range made of two integers separated from a comma 
#                               // remember first residue is number 1  do not put zero
#params.special_helix = false   // this parameter is used in the phase of colouring for differential colouring, it can be true or false
#params.create_sh_ref = false   // used to transfer the annotation from a reference sequence onto others that are well aligned to it
*/

// retrival of fasta sequences on the basis of a one per line header-id file
// from n target files manually specified 

// #files_to_search = '/home/alessio/Desktop/erasmus-internship/backups/MIA2_vert.fasta' // '/home/alessio/Desktop/erasmus-internship/backups/all_mia3_got_so_far_oct14.fasta,/home/alessio/Desktop/erasmus-internship/backups/only_Vertebrata.fasta,/home/alessio/Desktop/erasmus-internship/backups/MIA3_from_drome_danre_homo.fasta'   to change    if needed


process nfiles_search {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${outpath}"}

	input:
	path inheader
	val files_to_search
	path py_script1

	output:
	path "${outpath}", emit: retrieved_fasta
	//stdout emit: standardout
	
	script:
	outpath = "${inheader}".split('\\.')[0] + ".fasta"
	"""
	python3 ${py_script1} ${inheader} ${files_to_search} ${outpath}
	"""
}



// In this section if the --TRIMM option is given from command line the sequences in input
// will be cut right untill n residues before the last residue with posterior prob greater than 0.9
// where n is the argument of the trim option


process trimm_seqs {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
        container params.CONTAINER
        tag { "${trimmed_file}"}

	input:
	path in_multifasta
	path py_script2
	val trimm_val

	output:
	path "${trimmed_file}", emit: trimmed_fasta
	//stdout emit: standardout
	

	script:
	trimmed_file = "${in_multifasta}".split('\\.')[0] + "-trimm.fasta"
	"""
	python3 ${py_script2} ${in_multifasta} ${params.TRIMM_PLP_DIR} ${trimmed_file} ${trimm_val} ${params.SIGNALPEPT}
	"""
}



// In this section an allignment is generated on the basis of the fassta obtained
// in the previous process, the html file is deleted right after
// because a better one will be generated again later 

process align_generation {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container "cbcrg/tcoffee@sha256:36c526915a898d5c15ede89bbc3854c0a66cef22db86285c244b89cad40fb855" 
	tag { "${tcoffee_outfilepath}" }
	afterScript "rm ${tcoffee_otfilepath}.html"

	input:
	path multifasta
	
	output:
        path tcoffee_outfilepath, emit: aln_file

	script:
	tcoffee_outfilepath = "${multifasta}".split('\\.')[0] + '.aln'
	"""
	export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/ncbi-blast/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/ TEMP=/tmp PERL5LIB=/opt/tcoffee/perl/lib/perl5 DIR_4_TCOFFEE=/opt/tcoffee EMAIL_4_TCOFFEE=tcoffee.msa@gmail.com CACHE_4_TCOFFEE=/tmp/cache/ LOCKDIR_4_TCOFFEE=/tmp/lck/ TMP_4_TCOFFEE=/tmp/				# for singularity env variables
	t_coffee -in ${multifasta} -outfile ${tcoffee_outfilepath}
	"""
}



// In this process the rename file needed from t-coffee is created using the alignment file just created
// the file will be like, name present in the aln file for that sequence, space, new name for that sequence 
// the length of the new names is set to be 15 that's maximum length for sequence names in t-coffee html color 
// and the new name are the ids and scientific names 

process rename_gen {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${outfile_path}" }
	//afterScript "rm -f ${color_list_to_remove}"    this is a trick to force nextflow to recreate the color file each time otherwise it will just reuse it without recreating it since it sees no change to the last time modified this behaviour can lead to tricky problems

	input:
	path aln_outfile 
	
	output:
	path outfile_path, emit: rename_file
	stdout emit: standardout     

	script:
	outfile_path = "${aln_outfile}".split('\\.')[0] + '.rename'
	//color_list_to_remove = "${aln_outfile}".split('\\.')[0] + '.color'   this has indeed the same name of what will be created after so it is wiped out every time right before the step of creation of such file
	"""
	#!/usr/bin/env python3
	
	#print('$aln_outfile')
	
	with open('$aln_outfile', 'r') as in_aln, open('$outfile_path', 'w') as out_rename:
		title = in_aln.readline()
		block = -1
		while block < 1:
			line = in_aln.readline()
			if line[0] == '\\n' or line[0] == ' ':
				block += 1
			else:
				old_name = line.split(' ')[0]
				new_name = ''
				if '1_1' in old_name:
					splitted = old_name.split('_')
					#new_name = splitted[0] + '-' + splitted[1] +  '_'  + old_name.split('1_1_')[1].split('-')[0] + '                       '
					new_name = old_name.split('-')[-1] + '.'  + old_name.split('1_1_')[1].split('-')[0] + '                       '
				elif '!' in old_name:
					#new_name = old_name.split('-')[3] +'-'+ old_name.split('-')[4] + '_' + (old_name.split('_')[0]).split('!')[1] + '                    '
					new_name = old_name.split('-')[-3] + '.' + (old_name.split('_')[0]).split('!')[1] + '                    '
				elif len(old_name.split('-')) >= 3:
					#new_name = old_name.split('-')[3] +'-'+ old_name.split('-')[4] +'_'+ old_name.split('-')[0] + '                  '
					new_name = old_name.split('-')[-3] +'.'+ old_name.split('-')[0] + '                  '
				else:
					new_name = old_name + '                  '
				actual_rule = old_name + ' ' + new_name[0:15] + '\\n'		# where the length of title is set
				#print(actual_rule)
				print(old_name)
				out_rename.write(actual_rule)
	"""
}



// this process basically creates a file that have on each line 
// aln old name, before the rename, <space> residue-position_in_sequence with posterior prob over threshold <space> color number from 0 to 9
// the colors are assigned on the basis of which range bin the pp value falls in, 
// the ranges are threshold_value 0,1 0,2 0,3 ecc


process residues_colors {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${aacolor_list}" }

	input: 
	each path(plp_file)
	path rename_out_file
	path color_aa_pyscript


	output:
	stdout emit: standardout		                             // used only for check during debugging
	path aacolor_list, emit: aacolor_list_outfile
	path specialH_aacolor_list, emit: specialH_aacolor_list_outfile     // used for the flag special_helix

	script:
	aacolor_list = "${plp_file}".split('\\.')[0].split('!')[-1] + '.color'
	specialH_aacolor_list = "${plp_file}".split('\\.')[0].split('!')[-1] + '.sHcolor'
	left_trimm = "${params.TRIMM}".split(',')[0]
	check2 = left_trimm.isNumber()
	if(params.TRIMM == false) {
		if(params.SPECIAL_HELIX == false) {
			"""
			python3 ${color_aa_pyscript} ${plp_file} ${rename_out_file} > ${aacolor_list}
			touch ${specialH_aacolor_list}
			"""
		} else {
			"""
			python3 ${color_aa_pyscript} ${plp_file} ${rename_out_file} > ${aacolor_list}
			python3 ${color_aa_pyscript} ${plp_file} ${rename_out_file} ${params.SPECIAL_HELIX} > ${specialH_aacolor_list}
			"""
		}
	} else if( "$check2"=="true" ) {
		if(params.SPECIAL_HELIX == false) {
			"""
			python3 ${color_aa_pyscript} ${plp_file} ${rename_out_file} ${params.TRIMM} ${params.SIGNALPEPT} > ${aacolor_list}
			touch ${specialH_aacolor_list}
			"""
		} else {
			"""
			python3 ${color_aa_pyscript} ${plp_file} ${rename_out_file} ${params.TRIMM} ${params.SIGNALPEPT} > ${aacolor_list}
			python3 ${color_aa_pyscript} ${plp_file} ${rename_out_file} ${params.TRIMM} ${params.SIGNALPEPT} ${params.SPECIAL_HELIX} > ${specialH_aacolor_list}
			"""
		}
	} else {
		"""
		echo '--TRIMM argument has been given but it does not seem to be an integer'   
		"""
	}
}



// Couloring the residues in the allignment based on the Hydrophobicity scale, default Kyte-Doolittle
// take a look at hydrophobicity_color_list_tcoffee_prepare.py  for the list of possible scales

process residues_hydrophobicity_colors {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
        container params.CONTAINER
        tag { "${hydrocolor_list_outfilepath}" }

	input:
	path in_fasta
	path rename_outfile 
	path hydro_color_aa_pyscript

	output:
	stdout emit: standardout            // used only for check during debugging
	path hydrocolor_list_outfilepath, emit: hydrocolor_list_file

	script:
	hydrocolor_list_outfilepath = "${rename_outfile}".split('\\.')[0] + "-${params.HYDRO_SCALE}.hydrocolor" 
	if(params.MASK == false)
		"""
		python3 ${hydro_color_aa_pyscript} ${in_fasta} ${rename_outfile} ${hydrocolor_list_outfilepath} ${params.HYDRO_SCALE}
		"""
	else
		"""
		python3 ${hydro_color_aa_pyscript} ${in_fasta} ${rename_outfile} ${hydrocolor_list_outfilepath} ${params.HYDRO_SCALE} ${params.MASK}
		"""
}


// The following processisalater update of the pipeline to handle also when it isfeda phobius -short prediction file 


process phob_mode_checker {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
        container params.CONTAINER
        tag { "phob_mode_checker" }

	input:
        path phob_stdout_file

	output:
	stdout emit: standardout

	script:
	"""
	if grep -q 'FT   DOMAIN' ${phob_stdout_file}; then
		echo 'long'
	else
		echo 'short'
	fi
	"""
}


// The following process creates another identical coloured alignment based on phobius stdoutput default format (-long).
// It simply Highlights the coloured segment predicted.


process phob_stout_colours {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container "alessiovignoli3/tango-project@sha256:36cc270916232308969735637dba81b775916b2d221a811ec13dec597296fe0b" // field retriever
	tag { "${phob_stout_colourlist}" }
	

	input:
	path plp_folder
        path rename_filepath 
	path phobius_stdout_filepath
	path pyscript
	val mode_phobius

	output:
	stdout emit: standardout            // used only for check during debugging
	path phob_stout_colourlist, emit: phob_stout_colourfile

	script:
	phob_stout_colourlist = "${rename_filepath}".split('\\.')[0] + '.tmcolor'
	mode_phob = "${mode_phobius}".trim()
	if(params.TRIMM == false)
		if(params.SPECIAL_HELIX == false)
			"""
			python3 ${pyscript} ${phobius_stdout_filepath} ${phob_stout_colourlist} ${rename_filepath} ${mode_phob}
			"""
		else
			"""
			python3 ${pyscript} ${phobius_stdout_filepath} ${phob_stout_colourlist} ${rename_filepath} ${mode_phob} ${params.SPECIAL_HELIX} 
			"""
	else
		if(params.SPECIAL_HELIX == false)
			"""
			python3 ${pyscript} ${phobius_stdout_filepath} ${phob_stout_colourlist} ${rename_filepath} ${mode_phob} ${plp_folder} ${params.TRIMM} ${params.SIGNALPEPT}
			"""
		else
        	        """
			python3 ${pyscript} ${phobius_stdout_filepath} ${phob_stout_colourlist} ${rename_filepath} ${mode_phob} ${plp_folder} ${params.TRIMM} ${params.SIGNALPEPT} ${params.SPECIAL_HELIX}
			"""
}


// In this step an html coloured file is generated
//  for this reason this process depends on all the four above


process couloring_aln  {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false, saveAs: { filename -> if (params.SPECIAL_HELIX != false) filename
										else if (filename.endsWith("-specialH_pp.html")) null
										else filename
										})
	container "cbcrg/tcoffee@sha256:36c526915a898d5c15ede89bbc3854c0a66cef22db86285c244b89cad40fb855" //Version_13.45.47.aba98c5
	tag { "${html_outfile}" }
	
	input:
	path aln_file
	path rename_file 
	path color_file
	path sHcolor_file
	path hydrocolor_file 
	path phob_tmcolor_file 

	output:
	stdout emit: standardout
	path "*.html", emit: coulor_out_files
	path ascii_scorefile, emit: ascii_score
	
	script:
	html_outfile = "${rename_file}".split('\\.')[0] + '-pp.html'
	html_sH_outfile = "${rename_file}".split('\\.')[0] + '-specialH_pp.html'
	renamed_aln = "${rename_file}".split('\\.')[0] + '-renamed.aln'
	html_hydro_outfile = "${rename_file}".split('\\.')[0] + "-${params.HYDRO_SCALE}-hydro.html"
	html_phob_tmcolor_outfile = "${rename_file}".split('\\.')[0] + '-phobtm.html'
	ascii_scorefile = "${rename_file}".split('\\.')[0] + '.score_ascii'
	if(params.SPECIAL_HELIX == false)
		"""
		export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/ncbi-blast/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/ TEMP=/tmp PERL5LIB=/opt/tcoffee/perl/lib/perl5 DIR_4_TCOFFEE=/opt/tcoffee EMAIL_4_TCOFFEE=tcoffee.msa@gmail.com CACHE_4_TCOFFEE=/tmp/cache/ LOCKDIR_4_TCOFFEE=/tmp/lck/ TMP_4_TCOFFEE=/tmp/                           # for singularity env variables
		export ALN_LINE_LENGTH=150		# to make allignments of that length in html format
		touch ${html_sH_outfile}
		echo ${html_outfile}
		t_coffee -other_pg seq_reformat -in ${aln_file} -rename ${rename_file} -out ${renamed_aln}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +evaluate blosum62mt -output score_ascii > ${ascii_scorefile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${color_file} -output color_html > ${html_outfile}
		sed -i 's/9B92FF/FFFF99/' ${html_outfile}    # all of this do change every line that have 
		sed -i 's/B4FFB4/FFFF00/' ${html_outfile}    # the pattern specified
		sed -i 's/BEFFBE/CDFF00/' ${html_outfile}
		sed -i 's/C8FFC8/00FF3C/' ${html_outfile}
		sed -i 's/FFFFB7/00FF91/' ${html_outfile}
		sed -i -z 's/FFFFAD/00FFBC/' ${html_outfile} # this one changes only the first occurrence
		sed -i 's/FFFFAD/00FFF7/' ${html_outfile}    # note this has the same 1 field as above
		sed -i 's/FFE6E6/00E6FF/' ${html_outfile}
		sed -i 's/FFDCDC/00C4FF/' ${html_outfile}
		sed -i 's/FFC8C8/00ABFF/' ${html_outfile}
		echo ${html_hydro_outfile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${hydrocolor_file} -output color_html > ${html_hydro_outfile}
		sed -i 's/9B92FF/6666FF/' ${html_hydro_outfile} # changing color scheme of hydrophob file
		sed -i 's/B4FFB4/66B2FF/' ${html_hydro_outfile}
		sed -i 's/BEFFBE/99CCFF/' ${html_hydro_outfile}
		sed -i 's/C8FFC8/CCE5FF/' ${html_hydro_outfile}
		sed -i 's/FFFFB7/99FFFF/' ${html_hydro_outfile}
		sed -i -z 's/FFFFAD/FFCCFF/' ${html_hydro_outfile}
		sed -i 's/FFFFAD/FFE5CC/' ${html_hydro_outfile}
		sed -i 's/FFE6E6/FFCCCC/' ${html_hydro_outfile}
		sed -i 's/FFDCDC/FF9999/' ${html_hydro_outfile}
		sed -i 's/FFC8C8/FF6666/' ${html_hydro_outfile}
		echo ${html_phob_tmcolor_outfile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${phob_tmcolor_file} -output color_html > ${html_phob_tmcolor_outfile}
		rm -f ${renamed_aln} ${color_file} ${hydrocolor_file} ${phob_tmcolor_file} #Removing all the non output files. the files necessary for the scripts so far but not for an output point of view
		"""
	else
		"""
		export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/ncbi-blast/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/ TEMP=/tmp PERL5LIB=/opt/tcoffee/perl/lib/perl5 DIR_4_TCOFFEE=/opt/tcoffee EMAIL_4_TCOFFEE=tcoffee.msa@gmail.com CACHE_4_TCOFFEE=/tmp/cache/ LOCKDIR_4_TCOFFEE=/tmp/lck/ TMP_4_TCOFFEE=/tmp/                           # for singularity env variables
		export ALN_LINE_LENGTH=150              # to make allignments of that length in html format
		echo ${html_outfile}
		t_coffee -other_pg seq_reformat -in ${aln_file} -rename ${rename_file} -out ${renamed_aln}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +evaluate blosum62mt -output score_ascii > ${ascii_scorefile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${color_file} -output color_html > ${html_outfile}
		sed -i 's/9B92FF/FFFF99/' ${html_outfile}    # all of this do change every line that have 
		sed -i 's/B4FFB4/FFFF00/' ${html_outfile}    # the pattern specified
		sed -i 's/BEFFBE/CDFF00/' ${html_outfile}
		sed -i 's/C8FFC8/00FF3C/' ${html_outfile}
		sed -i 's/FFFFB7/00FF91/' ${html_outfile}
		sed -i -z 's/FFFFAD/00FFBC/' ${html_outfile} # this one changes only the first occurrence
		sed -i 's/FFFFAD/00FFF7/' ${html_outfile}    # note this has the same 1 field as above
		sed -i 's/FFE6E6/00E6FF/' ${html_outfile}
		sed -i 's/FFDCDC/00C4FF/' ${html_outfile}
		sed -i 's/FFC8C8/00ABFF/' ${html_outfile}
		echo ${html_sH_outfile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${sHcolor_file} -output color_html > ${html_sH_outfile}
		sed -i 's/9B92FF/FFFF99/' ${html_sH_outfile}    # all of this do change every line that have 
		sed -i 's/B4FFB4/FFFF00/' ${html_sH_outfile}    # the pattern specified
		sed -i 's/BEFFBE/CDFF00/' ${html_sH_outfile}
		sed -i 's/C8FFC8/00FF3C/' ${html_sH_outfile}
		sed -i 's/FFFFB7/00FF91/' ${html_sH_outfile}
		sed -i -z 's/FFFFAD/00FFBC/' ${html_sH_outfile} # this one changes only the first occurrence
		sed -i 's/FFFFAD/00FFF7/' ${html_sH_outfile}    # note this has the same 1 field as above
		sed -i 's/FFE6E6/00E6FF/' ${html_sH_outfile}
		sed -i 's/FFDCDC/00C4FF/' ${html_sH_outfile}
		sed -i 's/FFC8C8/00ABFF/' ${html_sH_outfile}
		echo ${html_hydro_outfile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${hydrocolor_file} -output color_html > ${html_hydro_outfile}
		sed -i 's/9B92FF/6666FF/' ${html_hydro_outfile} # changing color scheme of hydrophob file
		sed -i 's/B4FFB4/66B2FF/' ${html_hydro_outfile}
		sed -i 's/BEFFBE/99CCFF/' ${html_hydro_outfile}
		sed -i 's/C8FFC8/CCE5FF/' ${html_hydro_outfile}
		sed -i 's/FFFFB7/99FFFF/' ${html_hydro_outfile}
		sed -i -z 's/FFFFAD/FFCCFF/' ${html_hydro_outfile}
		sed -i 's/FFFFAD/FFE5CC/' ${html_hydro_outfile}
		sed -i 's/FFE6E6/FFCCCC/' ${html_hydro_outfile}
		sed -i 's/FFDCDC/FF9999/' ${html_hydro_outfile}
		sed -i 's/FFC8C8/FF6666/' ${html_hydro_outfile}
		echo ${html_phob_tmcolor_outfile}
		t_coffee -other_pg seq_reformat -in ${renamed_aln} -action +color_residue ${phob_tmcolor_file} -output color_html > ${html_phob_tmcolor_outfile}
		rm -f ${renamed_aln} ${color_file} ${sHcolor_file} ${hydrocolor_file} ${phob_tmcolor_file} #Removing all the non output files. the files necessary for the scripts so far but not for an output point of view
		"""
}


/*

// Special section were the reference sequence annotation is transferred to other sequences. The reference sequence is the last one from the top in the allignment. The annotation transferred is the presence of phobius standard output prediction, meaning that the aa aligned with the start and end of the region of interest in the reference will be found , and along with the aa also the position of sauch aa on the sequence will be found.
// Basically this process founds where is the start and end position of a specific region along each sequences that do not have such prediction but are still aligned to the reference. 

tranfer_sh_annotation_pyscript =  scripts + "transfer_annotation.py"

process tranfer_sh_annotation {

	input:
	val aln_f from couloring_aln_out2
	val plp_folder
	val phobius_stdout_filepath
	path tranfer_sh_annotation_pyscript

	output:
	stdout tranfer_sh_annotation_out1

	shell:
	transfered_phobtm_pred = "$aln_f".split('\\.')[0] + '-' + params.trimm + '.transf_phobtm'
	if(params.create_sh_ref == false)
		"""
		echo 'the option --create_sh_ref has not been given no tranfer of annotation will be done'
		rm -f !{aln_f}
		"""
	else if(params.trimm == false)
		"""
		python3 !{tranfer_sh_annotation_pyscript} !{aln_f} !{phobius_stdout_filepath} > $transfered_phobtm_pred
		rm -f !{aln_f}
		"""
	else
		"""
		python3 !{tranfer_sh_annotation_pyscript} !{aln_f} !{phobius_stdout_filepath} !{plp_folder} !{params.trimm} !{params.signalpept} > $transfered_phobtm_pred
		rm -f !{aln_f}
		"""
}
tranfer_sh_annotation_out1.view()
*/




workflow retrieve_fastas {

	take:
	pattern_to_headers
	pattern_to_searchfiles

	main:
	in_headers = Channel.fromPath(pattern_to_headers)
	retrieval_pyscript = params.SCRIPTS + 'from_header_to_fasta.py'
	nfiles_search(in_headers, pattern_to_searchfiles, retrieval_pyscript)
	suffix = "fastaoneline"
	retrieved_file = oneliner_ch(nfiles_search.out.retrieved_fasta, suffix)

	emit:
	//stout = nfiles_search.out.standardout
	retrieved_file
}

workflow alignment_generation { 

	take:
	pattern_to_input
	
	main:
	in_aln = pattern_to_input
	prefix = "${params.INPUT}".split('/')[-1].split('\\.')[0]
	if ( params.TRIMM != false ) {
		trimm_pyscript = params.SCRIPTS + "trimm_multifasta.py"
		trimm_seqs(pattern_to_input, trimm_pyscript, params.TRIMM)
		in_aln = trimm_seqs.out.trimmed_fasta
		prefix = prefix + "-trimm"
	}
	align_generation(in_aln)
	rename_gen(align_generation.out.aln_file)
	Channel.fromPath(params.PLP_DIR, type:'dir').set{ plp_dir  }
	rename_gen.out.standardout.splitText().flatMap{ "${it}".trim() + ".plp"  }.set{ fasta_old_headers }
	plp_dir.combine(fasta_old_headers).flatMap{ "${it[0]}" + '/' + "${it[1]}" }.set{ plp_filepath }
	color_aa_pyscript = params.SCRIPTS + "pp_color_list_tcoffee_prepare.py"
	residues_colors(plp_filepath, rename_gen.out.rename_file, color_aa_pyscript)
	converger(residues_colors.out.aacolor_list_outfile, prefix, "color")
	
	emit:
	stout = converger.out.stout						//for debug
	sh_color_list = residues_colors.out.specialH_aacolor_list_outfile
	color_file = converger.out.converged_file
	prefix
	in_aln
	rename_file = rename_gen.out.rename_file
	plp_dir
	alignment = align_generation.out.aln_file
}

workflow sh_colouring {

	take:
	color_files
	prefix_name

	main:
	converger(color_files, prefix_name, "sHcolor")

	emit:
	shcolor_file = converger.out.converged_file
	stout = converger.out.stout
	
}

workflow allignment_couloring {

	take:
	in_fastas
	rename_out
	plp_folder
	aln_outfile
	pp_color_file
	pp_shcolor_file

	main:
	hydro_color_aa_pyscript = params.SCRIPTS + "hydrophobicity_color_list_tcoffee_prepare.py"
	residues_hydrophobicity_colors(in_fastas, rename_out, hydro_color_aa_pyscript)
	phob_stout_colour_pyscript = params.SCRIPTS + "phobiuos_stdout_color_list_tcoffee_prepare.py"
	Channel.fromPath(params.PHOB_STDOUT).set{phob_stdout}
	phob_mode_checker(phob_stdout)
	phob_stout_colours(plp_folder, rename_out, phob_stdout, phob_stout_colour_pyscript, phob_mode_checker.out.standardout)
	couloring_aln(aln_outfile, rename_out, pp_color_file, pp_shcolor_file, residues_hydrophobicity_colors.out.hydrocolor_list_file, phob_stout_colours.out.phob_stout_colourfile)

	emit:
	stout = couloring_aln.out.standardout
}

workflow {
	aln_gen = false
	if ( params.SEARCHFILES != false ) {
		retrieve_fastas(params.INPUT, params.SEARCHFILES)
		aln_gen = alignment_generation(retrieve_fastas.out.retrieved_file)
	} else {
		in_file = Channel.fromPath(params.INPUT)
		aln_gen = alignment_generation(in_file)
	}
	//retrieve_fastas.out.stout.view()
	//retrieve_fastas.out.retrieved_file.view()
	//aln_gen.stout.view()
	specialh_colour = Channel.fromPath(params.PHOB_STDOUT)		// it has to point to a file that surely exist in case below is false
	if ( params.SPECIAL_HELIX != false ) {
		sh_colouring(aln_gen.sh_color_list, aln_gen.prefix)
		//sh_colouring.out.stout.view()
		specialh_colour = sh_colouring.out.shcolor_file
	}
	allignment_couloring(aln_gen.in_aln, aln_gen.rename_file, aln_gen.plp_dir, aln_gen.alignment, aln_gen.color_file, specialh_colour)
	allignment_couloring.out.stout.view()
}

