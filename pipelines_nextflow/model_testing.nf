#!/usr/bin/env nextflow


// ### PLEASE CHECK THE COMMENTS AND COMMENTED LINES IN THE SCRIPTS or SHELLS BLOCKS ###


// tm = tranmembrane
// om = Original model
// M3 = mark three third model created of the phobius model 
// ns = negatie set
// ps = positive set TANGO1 proteins not used for the training of the model
// fs = frequency setthe TANGO1 proteins used for the assestment of the frequency of aa this prot are excluded from the positive set
// pred = prediction


// Where all python scripts are for this project

scripts = "/home/alessio/Desktop/erasmus-internship/scripts/"

// Test folder where the actual pipeline is launched most of the times and where some file along this pipeline are saved.

test_folder = "/home/alessio/Desktop/erasmus-internship/test/"

// The folder containning all the results obtained from phobious. having subdirectories named after the specific set and model used for the prediction (like test1). This subdirectories have each two subdirectory named always pngs and plps. The pngs subdir contains the pngs of the sequence in the set of the named parent dyrecory (test1), same for the plps subdir. This plp and png files have to be computed before this script. Another file that has to be present before the start of the script is the stdoutput of phobius for each of the three sets, again this file has to be in the respective directory (the parent dir of the pngs,plps that is the same as the subdir of the folder containning all the results obtained from phobious) (test1). This file can be created simply re-directing the stdout of phobius into a file. The variable for the name of such file can be found iin the section of the process SOV_computation.

phobiuos_results_folder = "/home/alessio/Desktop/erasmus-internship/all_phobius_results/" //to change if needed remeber to put the slash at the end, this folder is parent of the below
phobiuos_om_tm_ns_result_folder = phobiuos_results_folder + "tm_negative_set1_original_model/"
phobiuos_om_nontm_ns_result_folder = phobiuos_results_folder + "nontm_negative_set2_original_model/"
phobiuos_om_ps_result_folder = phobiuos_results_folder + "MIA2_vertebrate_orthologs_OM/"  // positive set predictions are searched here
// all this above directories have to be created and filled of png and plp before the run of the script in the architecture specified above, rmember each of theese three dir have pngs and plps subdir

// All multifasta and respective headers. Have to be present before the the lounch of the pipeline.

prefix = "_set" // the word preceding the set identifier, since there could be different ps,ns set trials each one is numbered
hext = ".headers"  // the extention of every header file
fext = ".fasta"

tm_ns_name = "tm_negative"           // used later on
nontm_ns_name = "nontm_negative"     // same
ps_name = "positive"                 // for
fs_name = "frequency_trainig"        // this

tm_ns_headers = test_folder + tm_ns_name + prefix + "1" + hext
tm_ns_multifasta = test_folder + tm_ns_name + prefix + "1" + fext
nontm_ns_headers = test_folder + nontm_ns_name + prefix + "2" + hext
nontm_ns_multifasta = test_folder + nontm_ns_name + prefix + "2" + fext
ps_headers = test_folder + ps_name + prefix + "3" + hext
ps_multifasta = test_folder + ps_name + prefix + "3" + fext
fs_headers = test_folder + fs_name + prefix + "3" + hext
fs_multifasta = test_folder + fs_name + prefix + "3" + fext

complete_suffix = "set1233_" // this will be used later and each number is composed by previously defined set numbers (lines above), could be 1211 if nontm_ns is set2 and the rest set 1

// Phobius directory

phobius_dir = "/home/alessio/phobius101_linux/tmp/tmpbyom7j/phobius/"

// Backups directory. Where the various models will be stored

backup_dir = "/home/alessio/Desktop/erasmus-internship/backups/"

// Parameters used in this script

params.mark = 'bubba'




// First process, copies the phobius-Mnumber specified by -mark option to the backups dir, if there is no file with such name in the test folder an error occurs. If such file exists it is superimposed to phobius.model file. Creating the actual new model

model_name = test_folder + "phobius-M" + params.mark
model_filename_usedby_phobius = phobius_dir + "phobius.model"

process check_and_backup {

	input:
	path model_name
	path backup_dir
	path model_filename_usedby_phobius
	
	output:
	stdout check_and_backup_out

	shell:
	"""
	echo 'This is the name of the model used:' !{model_name} 
	echo 'if an error like this -> montage-im6.q16: width or height exceeds limit <- occurs make shure to uncomment the code lines in the visual_check_gen section'
	cp !{model_name} !{backup_dir}
	cat !{model_name} > !{model_filename_usedby_phobius}
	"""
}
check_and_backup_out.view()



// In the below process the folder containing the outputs of the pipeline are created and named after thw names specified in the fasta,header section. Here the complete_suffix identifier is used to generate a unique folder 

process creation_of_dirs {

	input:
	val phobiuos_results_folder
	val complete_suffix
	val tm_ns_name
	val nontm_ns_name
	val ps_name

	output:
	val paths_list into creation_of_dirs_out1
	stdout creation_of_dirs_out2 
	val current_model_pred_dir into creation_of_dirs_out3

	script:
	current_model_pred_dir = phobiuos_results_folder + complete_suffix + "M" + params.mark + "/"
	subdir_tm_ns_name = current_model_pred_dir + tm_ns_name
	subdir_nontm_ns = current_model_pred_dir + nontm_ns_name
	subdir_ps_name = current_model_pred_dir + ps_name
	paths_list = current_model_pred_dir + " " + subdir_tm_ns_name + "/ " + subdir_nontm_ns + "/ " + subdir_ps_name + "/"
	"""
	mkdir -p $current_model_pred_dir
	mkdir -p $subdir_tm_ns_name
	mkdir -p $subdir_nontm_ns
	mkdir -p $subdir_ps_name
	"""
}
//creation_of_dirs_out1.view()




// In this process the prediction with the specified model is carried out. The directories created in the previous step are filled with posterior-prob png files. In adddition the stdout output of phobius that contains the actual prediction for every sequence is copied to a file. This output file will be then used for the computation of SOV measure.

process pred_step {

	input:
	val dirs_paths from creation_of_dirs_out1
	path tm_ns_multifasta
	path nontm_ns_multifasta
	path ps_multifasta

	output:
	stdout pred_step_out1
	val list_of_dirpaths into pred_step_out2

	script:
	list_of_dirpaths = dirs_paths.split(' ')
	parent_dir = list_of_dirpaths[0]
	tm_ns_dir = list_of_dirpaths[1]
	nontm_ns_dir = list_of_dirpaths[2]
	ps_dir = list_of_dirpaths[3]
	"""
	limit=`cat $tm_ns_multifasta | wc -l`
	i=2; while [ \$i -le \$limit ]; do
		head -n \$i $tm_ns_multifasta | tail -n 2 | phobius-M7.pl -png $tm_ns_dir`head -n \$i $tm_ns_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.png -plp $tm_ns_dir`head -n \$i $tm_ns_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.plp 
	i=\$((\$i + 2))
	done > $parent_dir'phobius_output_tm_ns.txt'

	limit=`cat $nontm_ns_multifasta | wc -l`
        i=2; while [ \$i -le \$limit ]; do
                head -n \$i $nontm_ns_multifasta | tail -n 2 | phobius-M7.pl -png $nontm_ns_dir`head -n \$i $nontm_ns_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.png -plp $nontm_ns_dir`head -n \$i $nontm_ns_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.plp
        i=\$((\$i + 2))
        done > $parent_dir'phobius_output_nontm_ns.txt'

	limit=`cat $ps_multifasta | wc -l`
        i=2; while [ \$i -le \$limit ]; do
                head -n \$i $ps_multifasta | tail -n 2 | phobius-M7.pl -png $ps_dir`head -n \$i $ps_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.png -plp $ps_dir`head -n \$i $ps_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.plp
        i=\$((\$i + 2))
        done > $parent_dir'phobius_output_ps.txt'
	"""
}
//pred_step_out1.view()
//pred_step_out2.view()



// In this following process the newly created png files are concatenated with the pngs obtained using the original phobius model. There will be three comparison pdf files, one for each set defined so far. Each file will have the same structure, two png per row, where the left one is the original phobius png prediction for protA and the right one is the current model png prediction for the same protA. This will make easier visual check and comparison having on the left the comparison and on the right thr nrw prediction. This first process will concatenate the right png filepaths in one line, in the right order for the command montage that will actually create the pdfs.

process intermidiate_step {

	input:
	val phobiuos_om_tm_ns_result_folder
	val phobiuos_om_nontm_ns_result_folder
	val phobiuos_om_ps_result_folder
	val current_model_dirs from pred_step_out2
	
	output:
	stdout intermidiate_step_out1

	script:
	om_tm_ns_dir = phobiuos_om_tm_ns_result_folder + "pngs/"
	om_nontm_ns_dir = phobiuos_om_nontm_ns_result_folder + "pngs/"
	om_ps_dir = phobiuos_om_ps_result_folder + "pngs/"
	parent_dir = current_model_dirs[0]
        tm_ns_dir = current_model_dirs[1]
        nontm_ns_dir = current_model_dirs[2]
        ps_dir = current_model_dirs[3]
	"""
	#!/usr/bin/env python3
	
	import os
		
	for om, M in [('$om_tm_ns_dir', '$tm_ns_dir'), ('$om_nontm_ns_dir', '$nontm_ns_dir'), ('$om_ps_dir', '$ps_dir')]:
		all_png_om = os.listdir(om)
		all_png_M = os.listdir(M)
		for png_filename in all_png_M:
			if png_filename.endswith(".png") and png_filename in all_png_om:
				if '!' in png_filename:
					correct_fiename = png_filename.split('!')[0] + '\\!' + png_filename.split('!')[1]
					print(om + correct_fiename + ' ' + M + correct_fiename + ' ', end='')
				else:
					print(om + png_filename + ' ' + M + png_filename + ' ', end='')
		print()
	"""

}
//intermidiate_step_out1.view()



// In this following process the actual pdf files used for visual check are created using the command montage that requires the paths to the files in a specific order, this has been done by the above process.

process visual_check_gen {

	input:
	val  one_line_filepaths from intermidiate_step_out1
	val output_folder from creation_of_dirs_out3

	output:
	stdout visual_check_gen_out1

	shell:
	tm_ns_filepaths = one_line_filepaths.split('\n')[0]
	nontm_ns_filepaths = one_line_filepaths.split('\n')[1]
	ps_filepaths = one_line_filepaths.split('\n')[2]
	"""
	# This following lines are used to modify a read only config file 
	# otherwise montage can not create very large files 
	# if you want to edit the value option of each line
	#sudo sed -i 's/policy domain="resource" name="height" value="16KP"/policy domain="resource" name="height" value="128MB"/' /etc/ImageMagick-6/policy.xml
	#sudo sed -i 's/policy domain="resource" name="width" value="16KP"/policy domain="resource" name="width" value="128MB"/' /etc/ImageMagick-6/policy.xml
	#sudo sed -i 's/policy domain="resource" name="area" value="128MB"/policy domain="resource" name="area" value="256MB"/' /etc/ImageMagick-6/policy.xml
	# this last command is used when montage has problem converting to pdf
	#sudo sed -i 's/  <policy domain="coder" rights="none" pattern="PDF" \\/>/<!-- <policy domain="coder" rights="none" pattern="PDF" \\/> -->/' /etc/ImageMagick-6/policy.xml
	montage -tile 2x -geometry +5+5 !{tm_ns_filepaths} !{output_folder}visual_check_tm.pdf
	montage -tile 2x -geometry +5+5 !{nontm_ns_filepaths} !{output_folder}visual_check_nontm.pdf
	montage -tile 2x -geometry +5+5 !{ps_filepaths} !{output_folder}visual_check_ps.pdf
	"""
}
//visual_check_gen_out1.view()



// As mentioned at the top of this file, (in the section where the phobius results directories are mapped to variables), the stdoutput long-format(default) of phobius has to be re-diredted to a file and saved in the respective set folder. Since a SOV measure is computed for all tm_ns_set, nontm_ns_set and ps_set there is the need for a reference for each, that is in fact the stdoutput of the original phobius model.

om_phobius_stdout_tm_ns = phobiuos_om_tm_ns_result_folder + "phobius_output_om_tm_ns.txt"
om_phobius_stdout_nontm_ns = phobiuos_om_nontm_ns_result_folder + "phobius_output_om_nontm_ns.txt"
om_phobius_stdout_ps = phobiuos_om_ps_result_folder + "phobius_output_om_ps.txt"
om_phobius_stdout_ps_fs_fake_sHelix_ref = phobiuos_om_ps_result_folder + 'phobius_output_om_ps-trimm-70,30.transf_phobtm'
SOV_computation_pyscript = scripts + "SOV_from_phobius_stdout.py" 

process SOV_computation {

	input:
	path om_phobius_stdout_tm_ns 
	path om_phobius_stdout_nontm_ns
	path om_phobius_stdout_ps
	path SOV_computation_pyscript
	path om_phobius_stdout_ps_fs_fake_sHelix_ref
	val pred_step_out2
	
	output:
	stdout SOV_computation_out1

	script:
	parent_dir = pred_step_out2[0]
	"""
	echo 'TM negative set SOV:'
	echo 'Transmembrane region'
	./$SOV_computation_pyscript $om_phobius_stdout_tm_ns $parent_dir'phobius_output_tm_ns.txt' TRANSMEM TRANSMEM
	echo  'Domain region (everything that is not transmembrane)'
	./$SOV_computation_pyscript $om_phobius_stdout_tm_ns $parent_dir'phobius_output_tm_ns.txt' DOMAIN DOMAIN
	echo  'special helix region (there should be no sequences ideally)'
	./$SOV_computation_pyscript $om_phobius_stdout_tm_ns $parent_dir'phobius_output_tm_ns.txt' TRANSMEM SPECIAL true
	echo -e '\nnon-TM negative set SOV:'
	echo 'Transmembrane region'
	./$SOV_computation_pyscript $om_phobius_stdout_nontm_ns $parent_dir'phobius_output_nontm_ns.txt' TRANSMEM TRANSMEM
	echo 'Domain region'
	./$SOV_computation_pyscript $om_phobius_stdout_nontm_ns $parent_dir'phobius_output_nontm_ns.txt' DOMAIN DOMAIN
	echo  'special helix region (there should be no sequences ideally)'
        ./$SOV_computation_pyscript $om_phobius_stdout_nontm_ns $parent_dir'phobius_output_nontm_ns.txt' TRANSMEM SPECIAL true
	echo -e '\npositive set SOV:'
	echo 'Transmembrane region Normal helix (the one well predicted from Original Model)'
	./$SOV_computation_pyscript $om_phobius_stdout_ps $parent_dir'phobius_output_ps.txt' TRANSMEM NORMAL
	echo 'Transmembrane region Special helix'
        ./$SOV_computation_pyscript $om_phobius_stdout_ps $parent_dir'phobius_output_ps.txt' TRANSMEM SPECIAL 22 # very important this number has to be computed by the user and it is the number of sequences that actually have the speciaL helix already predicted by phobius original model
	echo 'Special helix not predicted in Original Model region'
	./$SOV_computation_pyscript $om_phobius_stdout_ps_fs_fake_sHelix_ref $parent_dir'phobius_output_ps.txt' TRANSMEM TRANSMEM
	echo 'Domain region'
	./$SOV_computation_pyscript $om_phobius_stdout_ps $parent_dir'phobius_output_ps.txt' DOMAIN DOMAIN
	"""
}
SOV_computation_out1.view()

