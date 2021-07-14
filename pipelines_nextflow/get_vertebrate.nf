#!/usr/bin/env nextflow

// Where all python scripts are for this project

scripts = "${params.parent_ofall_dir}scripts/"

// where the OMA folder containing all oma files is

OMA_folder = "${params.parent_ofall_dir}OMA/"

// where the Metaphor folder containing all met files is

Met_folder = "${params.parent_ofall_dir}Metaphor/"

// Test folder where the actual pipeline is launched most of the times
// and where some file along this pipeline are saved, like union files

test_folder = "${params.parent_ofall_dir}test/"

// The folder containning all the results obtained from phobious

phobiuos_results_folder = "${params.parent_ofall_dir}all_phobius_results/${params.phobius_results_subdir}"



// All the parameters used in this script 

params.kingdom = 'Vertebrata'







// Filter and polish of OMA _raw files

OMA_raws = Channel.fromPath(OMA_folder + "*_raw.fasta")
filter_oma_script = (scripts + "filter_OMA.py")

process OMA_output_filename_generator {
	container 'python'

	input:
	path oma_raw from OMA_raws

	output:
	file oma_raw into OMA_output_filename_generator_out1
	stdout into OMA_output_filename_generator_out2

	"""
	#!/usr/bin/env python3

	old = '$oma_raw'
	lista = old.split('_raw')
	if len(lista) == 2:
		new = lista[0] + '.fasta'
	else:
		print('This file has _raw twice into the name, it should just be one the suffix at the end:\\n', old)
	#print('here is the old file:', old)
	#print('lista:', lista)
	#print('here is the new file:', new)
	print(new, end='')
        """
}

process OMA_refinment {
	container 'python'
	publishDir  "$OMA_folder", mode: 'copy', overwrite: false 

	input:
	file raw_oma from OMA_output_filename_generator_out1
	val filename_filtered_oma from OMA_output_filename_generator_out2
	path py_script from filter_oma_script

	output:
        file "${filename_filtered_oma}" into raw_and_polished_OMA
        
	script:
        """
        ./${py_script} ${raw_oma} ${filename_filtered_oma}
	"""

}
//raw_and_polished_OMA.view()





// Filter and polish of Met _raw files

Met_raws = Channel.fromPath(Met_folder + "*_raw.fasta")
filter_met_script = (scripts + "filter_Metaphor.py")

process Met_output_filename_generator {
        container 'python'

        input:
        path met_raw from Met_raws

        output:
        file met_raw into Met_output_filename_generator_out1
        stdout into Met_output_filename_generator_out2

        """
        #!/usr/bin/env python3

        old = '$met_raw'
        lista = old.split('_raw')
        if len(lista) == 2:
                new = lista[0] + '.fasta'
        else:
                print('This file has _raw twice into the name, it should just be one the suffix at the end:\\n', old)
        #print('here is the old file:', old)
        #print('lista:', lista)
        #print('here is the new file:', new)
        print(new, end='')
        """
}

process Met_refinment {
	container 'python'
        publishDir  "$Met_folder", mode: 'copy', overwrite: false

	input:
	file raw_met from Met_output_filename_generator_out1
        val filename_filtered_met from Met_output_filename_generator_out2
	path py_script from filter_met_script

	output:
	file "${filename_filtered_met}" into raw_and_polished_Met

	script:
	"""
	./${py_script} ${raw_met} ${filename_filtered_met}
	"""

}
//raw_and_polished_Met.view()





// Union section, here all the files obtained from the filterers steps are united
// in a single file. The script used also removes the fasta sequences that have one or more <X>

all_refined_OMA = raw_and_polished_OMA.reduce { a, b -> ; return "$a"+","+"$b" }
all_refined_Met = raw_and_polished_Met.reduce { a, b -> ; return "$a"+","+"$b" }
union_script = (scripts + "unionist.py")

process union {
	container 'python'
	
	input:
	val one_line_oma_paths from all_refined_OMA
	val one_line_met_paths from all_refined_Met
	//val test_folder
	path py_script from union_script

	output:
	file "*.fasta" into union_output

	shell:
	union_input = one_line_oma_paths.trim()+','+one_line_met_paths.trim()
	"""
	./!{py_script} !{union_input}
	"""

}
//union_output.view()





// Filter out all the sequences that do not belong to <params.kingdom> this variable 
// Has to be set manually inputed from command line if needed to change, it is important that is exactly how
// ebi ena taxonomy cals the various classifications

https = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"
//final_output_filename = test_folder + "only_" + params.kingdom + ".fasta"

process species_name_retriever {
	container 'python'

	input:
	path union_file from union_output
	val https

	output:
        stdout  species_name_retriever_out

	"""
	#!/usr/bin/env python3
	
	with open('$union_file', 'r') as infile:
		#print(infile)
		string_https = ''
		for line in infile:
			if line[0] == '>':
				#print(line, end='')
				specie_name = (line.split('[')[1]).split(']')[0]
				#print(specie_name)
				http = ' ' + '$https' + specie_name.split(' ')[0] + '%20' + specie_name.split(' ')[1]
				if http in string_https:
					continue
				else:
					string_https +=	http
		print(string_https[1:], end='')
	"""

}
//species_name_retriever_out.view()

process ncbi_query {
	container 'alessiovignoli3/tango-project:bash_curl'

	input:
	val species_name_retriever_out

	output:
        file "species_informations.txt" into ncbi_query_out

	script:
	"""
	for i in ${species_name_retriever_out}
	do
		curl \$i >> species_informations.txt
	done

	"""
}
//ncbi_query_out.view()

process kingdom_filter {
	container 'python'

	input:
	path species_informations from ncbi_query_out
	//path species_informations from "/home/alessio/Desktop/erasmus-internship/test/species_informations.txt"
	val params.kingdom

	output:
	stdout kingdom_filter_out

	"""
	#!/usr/bin/env python3

	with open('$species_informations', 'r') as infile:
		#print(infile)
		taxid = ''
		speciename = ''
		for line in infile:
			if 'taxId' in line:
				taxid = (line.split('"')[3]).rstrip() #.split(',')[0]
				#print(taxid)
				#print(line, end='')
			elif 'scientificName' in line:
				speciename = (line.split('"')[3]).rstrip() #.split(',')[0]
				#print(speciename)
			elif 'lineage' in line and "$params.kingdom" in line:
				#print(line)
				print(speciename, taxid + ',', end='')
	"""
}
//kingdom_filter_out.view()

process final_step_filter {
	container 'python'
	publishDir  "$test_folder", mode: 'copy', overwrite: false

	input:
	path union_file from union_output
	//path final_output_filename
	//val path_to_only_kingdom from final_output_filename
	val dict_of_species from kingdom_filter_out

	output:
	file "*.fasta" into final_step_filter_out

	"""
	#!/usr/bin/env python3
	
	outfilename = "only_" + "$params.kingdom" + ".fasta"
	dict_of_species = {}
	list_species_taxid = '$dict_of_species'.split(',')
	list_species_taxid.pop((len(list_species_taxid)-1))
	for elem in list_species_taxid:
		specie = elem.split(' ')[0] + ' ' + elem.split(' ')[1]
		taxid = elem.split(' ')[2]
		dict_of_species[specie] = taxid
	#print(dict_of_species)
	with open('$union_file', 'r') as infile, open(outfilename, 'w') as outfile:
		#print(infile, outfile)
		status_check = 0
		for line in infile:
			if line[0] == '>':
				status_check = 0
				speciename = (line.split('[')[1]).split(']')[0]
				if speciename in dict_of_species:
					new_header = line.split('[')[0] + dict_of_species[speciename] +'-'+ (line.split('[')[1]).split(' ')[0] +'-'+ (line.split(']')[0]).split(' ')[1] +'\\n'
					#print(line, end='')
					outfile.write(new_header)
					status_check = 1
			elif status_check == 1:
				#print(line, end='')
				outfile.write(line)
	"""
}
//final_step_filter_out.view()



// Now for each sequence retrieved so far the postirior probability plot png
// and associated plp files will be generated using phobious program
// all the files will be saved in the <all_phobius_results> directory
// a momentary temporary will be created filled with one line fastas 
// for each speceis and than at the end removed 

png_dir = phobiuos_results_folder + "pngs/"
plp_dir = phobiuos_results_folder + "plps/"
tmp_dir = test_folder + "tmp/"


process all_in_one_phobius {
	container "alessiovignoli3/tango-project:phobius_image"
	beforeScript "mkdir -p ${phobiuos_results_folder} ${png_dir} ${plp_dir}"
	publishDir "${phobiuos_results_folder}", mode: 'move', overwrite: false, saveAs: { filename ->
			if (filename.endsWith(".png")) "pngs/$filename"
			else if (filename.endsWith(".plp")) "plps/$filename"
			else if (filename.endsWith(".txt")) filename
			else filename : null
		}

        input:
	path only_kingdom_multifasta from final_step_filter_out

        output:
        path "*.png" into all_in_one_phobius_out1
	path "*.plp" into all_in_one_phobius_out2
	path "*.txt" into all_in_one_phobius_out3

        script:
	"""
	export PATH="/home/phobius101_linux/tmp/tmpbyom7j/phobius:$PATH"
        limit=`cat $only_kingdom_multifasta | wc -l`
        i=2; while [ \$i -le \$limit ]; do
                head -n \$i $only_kingdom_multifasta | tail -n 2 | phobius.pl -png `head -n \$i $only_kingdom_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.png -plp `head -n \$i $only_kingdom_multifasta | tail -n 2 | head -n 1 | cut -d '>' -f 2`.plp
        i=\$((\$i + 2))
        done > 'phobius_output_only_${params.kingdom}.txt'
	"""

}




process one_specie_per_fasta {

	beforeScript "mkdir ${tmp_dir}"

	input:
	path input_fasta from final_step_filter_out

	output:
	stdout into one_specie_per_fasta_out

	"""
	#!/usr/bin/env python3
	
	#print('$tmp_dir')
	with open('$input_fasta', 'r') as infasta:
		#print(infasta)
		outfastaname = ''
		header = ''
		for line in infasta:
			if line[0] == '>':
				outfastaname = line.rstrip()[1:]
				header = line
			else:
				outfasta_filename = outfastaname
				path_to_outfasta = '$tmp_dir' + outfasta_filename + '.fasta'
				print(outfasta_filename, end=' ')
				with open(path_to_outfasta, 'w') as outfasta:
					outfasta.write(header)
					outfasta.write(line)
	
	"""
}
//one_specie_per_fasta_out.view()

process all_in_one_phobius {

	afterScript "rm -r ${tmp_dir}"

	input:
	val png_dir
	val plp_dir
	val tmp_dir
	val one_specie_per_fasta_out

	output:
	stdout all_in_one_phobius_out

	script:
	//len_line = len(one_specie_per_fasta_out)-1
	oneline_filepaths = one_specie_per_fasta_out.trim()
	"""
	for i in ${oneline_filepaths}
	do	
		phobius.pl -png ${png_dir}\$i.png -plp ${plp_dir}\$i.plp ${tmp_dir}\$i.fasta
	done

	"""
}
all_in_one_phobius_out.view()
*/
