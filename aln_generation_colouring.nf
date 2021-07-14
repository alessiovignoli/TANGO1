#!/usr/bin/env nextflow

// Where all python scripts are for this project

scripts = "${params.parent_ofall_dir}scripts/"

// Test folder where the actual pipeline is launched most of the times
// and where some file along this pipeline are saved, like union files

test_folder = "${params.parent_ofall_dir}test/"

// The folder containning all the results obtained from phobious for the sequences that we want to allign and colour. The sequences are coloured on the basis of the postirior probability file generated from phobius (.plp). This plp files have to be computed before this script.

phobiuos_results_folder = "${params.parent_ofall_dir}all_phobius_results/MIA2_vertebrate_orthologs_OM/"    // to change if needed remeber to put the slash at the end
phobius_stdout_filepath = phobiuos_results_folder + ""   // to change , this file is created from the simple re-direction of phobious stdout default, is neede in the phob_stout_colours section. The set originating this prediction has to be the set of sequences presnt in <input_header_file> input file

// The folder containing all plp files, posterior probability files generated from phobius

trimm_plp_folder ='/home/alessio/Desktop/erasmus-internship/all_phobius_results/good_pp_vert_originalModel/plps/'  // plp directory used in the trimm step
plp_folder = phobiuos_results_folder + "plps/"   // put here whatever the folder of the plp files is called

// Here there is the collection of the flags _ parameters specifiable from command line 
// each flag is called using two consecutive minus signs -- 
// and the name reported here after params 
// The signal peptide parameter is mandatory if --trimm option is specified

params.trimm = false           // this can be an integer or two intger separated by a comma whre the first one is cut on the left and the second on the right
params.signalpept = false      // this is an integers it has to go with trimm parameter
params.mask = false            // this is a range made of two integers separated from a comma 
                               // remember first residue is number 1  do not put zero
params.special_helix = false   // this parameter is used in the phase of colouring for differential colouring, it can be true or false
params.create_sh_ref = false   // used to transfer the annotation from a reference sequence onto others that are well aligned to it




// retrival of fasta sequences on the basis of a one per line header-id file
// from n target files manually specified 

input_header_file = test_folder + "MIA2_vert.headers"	// to change
files_to_search = '/home/alessio/Desktop/erasmus-internship/backups/MIA2_vert.fasta' // '/home/alessio/Desktop/erasmus-internship/backups/all_mia3_got_so_far_oct14.fasta,/home/alessio/Desktop/erasmus-internship/backups/only_Vertebrata.fasta,/home/alessio/Desktop/erasmus-internship/backups/MIA3_from_drome_danre_homo.fasta'   to change    if needed
output_path = test_folder + "MIA2_vert.fasta"    // to change
retrieval_pyscript = scripts + "from_header_to_fasta.py"

process nfiles_search {

	input:
	val input_header_file
	val files_to_search
	val output_path
	path retrieval_pyscript

	output:
	stdout nfiles_search_out

	shell:
	"""
	./!{retrieval_pyscript} !{input_header_file} !{files_to_search} !{output_path}
	"""
}
//nfiles_search_out.view()



// In this section if the -trimm option is given from command line the sequences in input
// will be cut right untill n residues before the last residue with posterior prob greater than 0.9
// where n is the argument of the trim option

trimm_pyscript = scripts + "trimm_multifasta.py"

process trimm_seqs {

	input:
	val nfiles_search_out
	path trimm_pyscript

	output:
	stdout trimm_seqs_out1
	val trimmed_file into trimm_seqs_out2

	script:
	in_multifasta = nfiles_search_out.trim()
	trimmed_file = in_multifasta.split('\\.')[0] + '-trimm.fasta'
	if(params.trimm == false)
		"""
		echo $in_multifasta
		"""
	else
		"""
		echo $trimmed_file
		./$trimm_pyscript $in_multifasta $trimm_plp_folder $trimmed_file $params.trimm $params.signalpept
		"""
}
//trimm_seqs_out.view()



// In this section an allignment is generated on the basis of the fassta obtained
// in the previous process, 6 cores are used and the html file is deleted right after
// because a better one will be generated again later 

process align_generation {

	afterScript "rm ${tcoffee_otfilepath}.html"

	input:
	val  trimm_seqs_out1
	
	output:
        val tcoffee_otfilepath into align_generation_out

	shell:
	tcoffee_infilepath = trimm_seqs_out1.trim()
	tcoffee_otfilepath = tcoffee_infilepath.split('\\.')[0] + '.aln'
	"""
	t_coffee -in !{tcoffee_infilepath} -n_core=6 -outfile !{tcoffee_otfilepath}
	"""
}
//align_generation_out.view()
				// when removing this also remove next commented script line



// In this process the rename file needed from t-coffee is created using the alignment file just created
// the file will be like, name present in the aln file for that sequence, space, new name for that sequence 
// the length of the new names is set to be 15 that's maximum length for sequence names in t-coffee html color 
// and the new name are the ids and scientific names 

process rename_gen {

	afterScript "rm -f ${color_list_to_remove}"   // this is a trick to force nextflow to recreate the color file each time otherwise it will just reuse it without recreating it since it sees no change to the last time modified this behaviour can lead to tricky problems

	input:
	val aln_infile from align_generation_out	//this one has to be removed to enable depemdencie from second process
	//val aln_infile from '/home/alessio/Desktop/erasmus-internship/test/frequency_trainig_set.aln'
	//val aln_infile from '/home/alessio/Desktop/erasmus-internship/test/frequency_trainig_set-trimm.aln'
	
	output:
	val outfile_path into rename_gen_out1
	stdout into rename_gen_out2                //used later on to search plp filenames in a broader dir

	script:
	outfile_path = aln_infile.split('\\.')[0] + '.rename'
	color_list_to_remove = aln_infile.split('\\.')[0] + '.color'  // this has indeed the same name of what will be created after so it is wiped out every time right before the step of creation of such file
	"""
	#!/usr/bin/env python3
	
	#print('$aln_infile')
	
	with open('$aln_infile', 'r') as in_aln, open('$outfile_path', 'w') as out_rename:
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
					new_name = splitted[0] + '-' + splitted[1] + '..............................'
				elif '!' in old_name:
					new_name = (old_name.split('_')[0]).split('!')[1] +'_'+ old_name.split('-')[3] +'-'+ old_name.split('-')[4] +'..............................'
				else:
					new_name = old_name.split('-')[0] +'_'+ old_name.split('-')[3] +'-'+ old_name.split('-')[4] +'..............................'
				actual_rule = old_name + ' ' + new_name[0:15] + '\\n'		# whre the length of title is set
				#print(actual_rule)
				print(old_name, end=',')
				out_rename.write(actual_rule)
	"""
}
//rename_gen_out1.view()
//rename_gen_out2.view()



// this process depends on the stdout output of the above and basically creates a file that have on each line 
// aln old name, before the rename, space residue-position_in_sequence with posterior prob over threshold space color number from 0 to 9
// the colors are assigned on the basis of which range bin the pp value falls in, 
// the ranges are threshold_value 0,1 0,2 0,3 ecc

plp_filenames = Channel.fromPath(plp_folder + "*.plp")
color_aa_pyscript = scripts + "pp_color_list_tcoffee_prepare.py"

process residues_colors {

	input:
	val plp_filepath from plp_filenames
	val rename_gen_out2
	val rename_out_file from rename_gen_out1
	path color_aa_pyscript

	output:
	stdout residues_colors_out1                             // used only for check during debugging
	val aacolor_list_outfilepath into residues_colors_out2
	val specialH_aacolor_list_outfilepath into residues_colors_out3    // used for the flag special_helix

	script:
	aacolor_list_outfilepath = rename_out_file.split('\\.')[0] + '.color'
	specialH_aacolor_list_outfilepath =rename_out_file.split('\\.')[0] + '.sHcolor'
	list_dirs = "$plp_filepath".split('/')
	len_list = list_dirs.size() - 1
	filename = list_dirs[len_list]
	len_non_extension = filename.size() - 5
	filename_without_ext = filename[0..len_non_extension]
	seq_id = filename_without_ext.split(':')[0]
	check = rename_gen_out2.contains(seq_id)
	left_trimm = "$params.trimm".split(',')[0]
	check2 = left_trimm.isNumber()
	if(params.trimm == false && "$check"=="true")
		if(params.special_helix == false)
			"""
			./$color_aa_pyscript $plp_filepath $rename_out_file >> $aacolor_list_outfilepath
			"""
		else
			"""
			./$color_aa_pyscript $plp_filepath $rename_out_file >> $aacolor_list_outfilepath
			./$color_aa_pyscript $plp_filepath $rename_out_file $params.special_helix >> $specialH_aacolor_list_outfilepath
			"""
	else if( "$check2"=="true" && "$check"=="true")
		if(params.special_helix == false)
			"""
			./$color_aa_pyscript $plp_filepath $rename_out_file $params.trimm $params.signalpept >> $aacolor_list_outfilepath
			"""
		else
			"""
			./$color_aa_pyscript $plp_filepath $rename_out_file $params.trimm $params.signalpept >> $aacolor_list_outfilepath
			./$color_aa_pyscript $plp_filepath $rename_out_file $params.trimm $params.signalpept $params.special_helix >> $specialH_aacolor_list_outfilepath
			"""

	else
		"""
		echo 'the file is not present'
		"""
}
//residues_colors_out1.view()
//residues_colors_out2.view()



// This process will create a color file based on a majority vote Hydrophobicity scale.
// this scale is deduced from a series of different scales such as Eisenberg and weisss, Engelman, Kyte Doolyttle, Hoop and Woods, Janin, Whimly White.
// this scale is arbitrary and deduced from the previously mentioned ones, the aa are coupled in pairs since there can be only 10 colors.
// the couples can be found inside the below script.

hydro_color_aa_pyscript = scripts + "hydrophobicity_color_list_tcoffee_prepare.py"

process residues_hydrophobicity_colors {

	input:
	path trimmed_fasta_filepath from trimm_seqs_out1
	val rename_out_filepath from rename_gen_out1
	path hydro_color_aa_pyscript

	output:
	stdout residues_hydrophobicity_colors_out1             // used only for check during debugging
	val hydrocolor_list_outfilepath into residues_hydrophobicity_colors_out2

	script:
	hydrocolor_list_outfilepath = rename_out_filepath.split('\\.')[0] + '.hydrocolor' 
	if(params.mask == false)
		"""
		./$hydro_color_aa_pyscript $trimmed_fasta_filepath $rename_out_filepath $hydrocolor_list_outfilepath
		"""
	else
		"""
		./$hydro_color_aa_pyscript $trimmed_fasta_filepath $rename_out_filepath $hydrocolor_list_outfilepath $params.mask
		"""
}
//residues_hydrophobicity_colors_out1.view()
//residues_hydrophobicity_colors_out2.view()



// The following process creates another identical coloured alignment based on phobius stdoutput default format (-long).
// It simply Highlights the coloured segment predicted.

phob_stout_colour_pyscript =  scripts + "phobiuos_stdout_color_list_tcoffee_prepare.py"

process phob_stout_colours {

	input:
	val plp_folder
        val rename_filepath from rename_gen_out1
	path phobius_stdout_filepath
	path phob_stout_colour_pyscript

	output:
	stdout phob_stout_colours_out1
	val phob_stout_colourlist into phob_stout_colours_out2

	script:
	phob_stout_colourlist = rename_filepath.split('\\.')[0] + '.tmcolor'
	if(params.trimm == false)
		if(params.special_helix == false)
			"""
			./$phob_stout_colour_pyscript $phobius_stdout_filepath $phob_stout_colourlist $rename_filepath
			"""
		else
			"""
			./$phob_stout_colour_pyscript $phobius_stdout_filepath $phob_stout_colourlist $rename_filepath $params.special_helix 
			"""
	else
		if(params.special_helix == false)
			"""
			./$phob_stout_colour_pyscript $phobius_stdout_filepath $phob_stout_colourlist $rename_filepath $plp_folder $params.trimm $params.signalpept
			"""
		else
        	        """
			./$phob_stout_colour_pyscript $phobius_stdout_filepath $phob_stout_colourlist $rename_filepath $plp_folder $params.trimm $params.signalpept $params.special_helix
			"""
}
phob_stout_colours_out1.view()



// In this step an html coloured file is generated
//  for this reason this process depends on all the four above

color_file = residues_colors_out2.last()
sHcolor_file = residues_colors_out3.last()
//aln_file = test_folder + "frequency_trainig_set.aln"		// this lines must be commented when the step of the 
//aln_file = test_folder + "frequency_trainig_set-trimm.aln"	// allignment is uncommented so as the dependencie below 

process couloring_aln  {
	
	input:
	val aln_file from align_generation_out     // this one is the dependencie cited before
	path rename_file from rename_gen_out1
	val rename_f from rename_gen_out1
	val color_file
	val sHcolor_file
	val hydrocolor_file from residues_hydrophobicity_colors_out2
	val phob_tmcolor_file from phob_stout_colours_out2
	val trimm_fasta_file from trimm_seqs_out2

	output:
	stdout couloring_aln_out1
	val aln_file into couloring_aln_out2

	shell:
	html_outfile = test_folder + "$rename_file".split('\\.')[0] + '-pp.html'
	html_sH_outfile = test_folder + "$rename_file".split('\\.')[0] + '-specialH_pp.html'
	renamed_aln = test_folder + "$rename_file".split('\\.')[0] + '-renamed.aln'
	html_hydro_outfile = test_folder + "$rename_file".split('\\.')[0] + '-hydro.html'
	html_phob_tmcolor_outfile = test_folder + "$rename_file".split('\\.')[0] + '-phobtm.html'
	ascii_scorefile = test_folder + "$rename_file".split('\\.')[0] + '.score_ascii'
	if(params.special_helix == false)
		"""
		echo !{html_outfile}
		t_coffee -other_pg seq_reformat -in !{aln_file} -rename !{rename_file} -out !{renamed_aln}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +evaluate blosum62mt -output score_ascii > !{ascii_scorefile}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{color_file} -output color_html > !{html_outfile}
		sed -i 's/9B92FF/FFFF99/' !{html_outfile}    # all of this do change every line that have 
		sed -i 's/B4FFB4/FFFF00/' !{html_outfile}    # the pattern specified
		sed -i 's/BEFFBE/CDFF00/' !{html_outfile}
		sed -i 's/C8FFC8/00FF3C/' !{html_outfile}
		sed -i 's/FFFFB7/00FF91/' !{html_outfile}
		sed -i -z 's/FFFFAD/00FFBC/' !{html_outfile} # this one changes only the first occurrence
		sed -i 's/FFFFAD/00FFF7/' !{html_outfile}    # note this has the same 1 field as above
		sed -i 's/FFE6E6/00E6FF/' !{html_outfile}
		sed -i 's/FFDCDC/00C4FF/' !{html_outfile}
		sed -i 's/FFC8C8/00ABFF/' !{html_outfile}
		echo !{html_hydro_outfile}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{hydrocolor_file} -output color_html > !{html_hydro_outfile}
		sed -i 's/9B92FF/6666FF/' !{html_hydro_outfile} # changing color scheme of hydrophob file
		sed -i 's/B4FFB4/66B2FF/' !{html_hydro_outfile}
		sed -i 's/BEFFBE/99CCFF/' !{html_hydro_outfile}
		sed -i 's/C8FFC8/CCE5FF/' !{html_hydro_outfile}
		sed -i 's/FFFFB7/99FFFF/' !{html_hydro_outfile}
		sed -i -z 's/FFFFAD/FFCCFF/' !{html_hydro_outfile}
		sed -i 's/FFFFAD/FFE5CC/' !{html_hydro_outfile}
		sed -i 's/FFE6E6/FFCCCC/' !{html_hydro_outfile}
		sed -i 's/FFDCDC/FF9999/' !{html_hydro_outfile}
		sed -i 's/FFC8C8/FF6666/' !{html_hydro_outfile}
		echo !{html_phob_tmcolor_outfile}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{phob_tmcolor_file} -output color_html > !{html_phob_tmcolor_outfile}
		rm -f !{renamed_aln} !{rename_f} !{color_file} !{hydrocolor_file} !{phob_tmcolor_file} !{trimm_fasta_file}  #Removing all the non output files. the files necessary for the scripts so far but not for an output point of view
		"""
	else
		"""
		echo !{html_outfile}
		t_coffee -other_pg seq_reformat -in !{aln_file} -rename !{rename_file} -out !{renamed_aln}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{color_file} -output color_html > !{html_outfile}
		sed -i 's/9B92FF/FFFF99/' !{html_outfile}    # all of this do change every line that have 
		sed -i 's/B4FFB4/FFFF00/' !{html_outfile}    # the pattern specified
		sed -i 's/BEFFBE/CDFF00/' !{html_outfile}
		sed -i 's/C8FFC8/00FF3C/' !{html_outfile}
		sed -i 's/FFFFB7/00FF91/' !{html_outfile}
		sed -i -z 's/FFFFAD/00FFBC/' !{html_outfile} # this one changes only the first occurrence
		sed -i 's/FFFFAD/00FFF7/' !{html_outfile}    # note this has the same 1 field as above
		sed -i 's/FFE6E6/00E6FF/' !{html_outfile}
		sed -i 's/FFDCDC/00C4FF/' !{html_outfile}
		sed -i 's/FFC8C8/00ABFF/' !{html_outfile}
		echo !{html_sH_outfile}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{sHcolor_file} -output color_html > !{html_sH_outfile}
		sed -i 's/9B92FF/FFFF99/' !{html_sH_outfile}    # all of this do change every line that have 
		sed -i 's/B4FFB4/FFFF00/' !{html_sH_outfile}    # the pattern specified
		sed -i 's/BEFFBE/CDFF00/' !{html_sH_outfile}
		sed -i 's/C8FFC8/00FF3C/' !{html_sH_outfile}
		sed -i 's/FFFFB7/00FF91/' !{html_sH_outfile}
		sed -i -z 's/FFFFAD/00FFBC/' !{html_sH_outfile} # this one changes only the first occurrence
		sed -i 's/FFFFAD/00FFF7/' !{html_sH_outfile}    # note this has the same 1 field as above
		sed -i 's/FFE6E6/00E6FF/' !{html_sH_outfile}
		sed -i 's/FFDCDC/00C4FF/' !{html_sH_outfile}
		sed -i 's/FFC8C8/00ABFF/' !{html_sH_outfile}
		echo !{html_hydro_outfile}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{hydrocolor_file} -output color_html > !{html_hydro_outfile}
		sed -i 's/9B92FF/6666FF/' !{html_hydro_outfile} # changing color scheme of hydrophob file
		sed -i 's/B4FFB4/66B2FF/' !{html_hydro_outfile}
		sed -i 's/BEFFBE/99CCFF/' !{html_hydro_outfile}
		sed -i 's/C8FFC8/CCE5FF/' !{html_hydro_outfile}
		sed -i 's/FFFFB7/99FFFF/' !{html_hydro_outfile}
		sed -i -z 's/FFFFAD/FFCCFF/' !{html_hydro_outfile}
		sed -i 's/FFFFAD/FFE5CC/' !{html_hydro_outfile}
		sed -i 's/FFE6E6/FFCCCC/' !{html_hydro_outfile}
		sed -i 's/FFDCDC/FF9999/' !{html_hydro_outfile}
		sed -i 's/FFC8C8/FF6666/' !{html_hydro_outfile}
		echo !{html_phob_tmcolor_outfile}
		t_coffee -other_pg seq_reformat -in !{renamed_aln} -action +color_residue !{phob_tmcolor_file} -output color_html > !{html_phob_tmcolor_outfile}
		rm -f !{renamed_aln} !{rename_f} !{color_file} !{sHcolor_file} !{hydrocolor_file} !{phob_tmcolor_file} !{trimm_fasta_file} #Removing all the non output files. the files necessary for the scripts so far but not for an output point of view
		"""
}
couloring_aln_out1.view()



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
		./!{tranfer_sh_annotation_pyscript} !{aln_f} !{phobius_stdout_filepath} > $transfered_phobtm_pred
		rm -f !{aln_f}
		"""
	else
		"""
		./!{tranfer_sh_annotation_pyscript} !{aln_f} !{phobius_stdout_filepath} !{plp_folder} !{params.trimm} !{params.signalpept} > $transfered_phobtm_pred
		rm -f !{aln_f}
		"""
}
tranfer_sh_annotation_out1.view()

