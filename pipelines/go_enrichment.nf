#!/usr/bin/env nextflow

// ### PLEASE CHECK THE COMMENTS AND COMMENTED LINES IN THE SCRIPTS or SHELLS BLOCKS ###


// Test folder where the actual pipeline is launched most of the times and where some file along this pipeline are saved.

test_folder = "${params.parent_ofall_dir}test/"




// Params used in this script

params.input_fastas_filenames = false      // defines whats the name or the glob path search to retrieve the input fasta no need 
					   // for .fasta extension is already present in the process
params.input_ids = false			// explained after
params.n_replicas = 3				// explained after


// The inputs of this pipeline are one or more fasta files containing the whole set on which the analysis has been done (whole_fasta_set),
// the ids (one per line) used in the headers of the fastas aformentioned outputted by the analysis done on which we do the GO enrichment search (input_ids)
// and the number of replicas of the GO enrichment analyses we want to do, with random selected ids of the input fastas, this will be the backgound to  which compare the results obtainedby GO on the input ids (params.n_replicas)

whole_fasta_set = Channel.fromPath(test_folder + "${params.input_fastas_filenames}.fasta")
input_ids = test_folder + "${params.input_ids}"



process ids_retrieval {
	container "python"

        input:
        path fastas from whole_fasta_set
	//path input_id_file from input_ids

        output:
        stdout ids_retrieval_out

        script:
        """
	#!/usr/bin/env python3
	
	with open('$fastas', 'r') as infasta:
		#print(infasta)
		for line in infasta:
			if '>' in line:
				header = (line.split(' ')[0]) 
				if '|' in header and 'OX' in line:
					uniprotid = header.split('|')[1]
					taxid = (line.split('OX')[1]).split(' ')[0][1:]
					print(uniprotid, taxid)
				elif '_' in header and 'TaxID' in line:
					uniprotid = header.split('_')[1]
					taxid = (line.split('TaxID')[1]).split(' ')[0][1:]
					print(uniprotid, taxid)
			else:
				continue
	"""
}
//ids_retrieval_out.view()
list_of_ids = ids_retrieval_out.toList()



process unifier {
	container 'ubuntu'
	//publishDir  "${test_folder}", mode: 'copy', overwrite: false

	input:
	val ids from list_of_ids

	output:
	//stdout into unifier_out
	file "all_uniprot_ids.txt" into unifier_out

	shell:	
	"""
	echo "!{ids}" > all_uniprot_ids.txt
	"""
}
//unifier_out.view()




def range = 1..params.n_replicas

process random_batch_creator {
	container 'ubuntu'
	//publishDir  "${test_folder}", mode: 'copy', overwrite: false

	input:
        path whole_ids from unifier_out
	path interesting_ids from input_ids
	val batch_number from range

        output:
        //stdout into random_batch_creator_out1
        file "random_batch_*.txt" into random_batch_creator_out1
	//file "!{interesting_ids}" into random_batch_creator_out2

	shell:
	"""
	batch_size=`wc -l !{interesting_ids} | cut -d ' ' -f 1`
	#echo \$batch_size
	sort -uR !{whole_ids} | head -n \$batch_size > random_batch_!{batch_number}.txt
	#| awk '{print}' ORS=' ' 
	#> random_batch_!{batch_number}.txt
	"""
}
random_batch_creator_out1.view()




process input_id_file_taxids_retrieval {
	container "python"
	publishDir  "${test_folder}", mode: 'copy', overwrite: false

	input:
	path total_ids from unifier_out
	path input_id_file from input_ids

	output:
	//stdout input_id_file_taxids_retrieval_out
	file "*.id_taxid" into input_id_file_taxids_retrieval_out

	script:
	"""
	#!/usr/bin/env python3

	output_file = "$input_id_file".split('.')[0] + ".id_taxid"
	with open('$total_ids', 'r') as totalset, open('$input_id_file', 'r') as infile, open(output_file, 'w') as outfile:
		#print(totalset, infile, outfile)
		one_line_headers = ''
		for in_line in infile:
			in_id = (in_line.rstrip()).split(' ')[0] + ' '
			one_line_headers += in_id
		#print(one_line_headers)
		for tot_line in totalset:
			tot_id = tot_line.split(' ')[0] + ' ' 
			if tot_id in one_line_headers:
				#print(tot_line, end='')
				outfile.write(tot_line)
        """
}
input_id_file_taxids_retrieval_out.view()



fai l exception block process che perogni id con 1 o 2 come taxid retriva i veri taxid magari usando il blocco qui sotto per l api
e poi integra tutto quello che esce da i processi ne species grouper che hai scritto venerdi 7 maggio

/*


process chunk_splitter {
	container 'python'

	input:
	path interest_ids from input_ids
        path batch_file from random_batch_creator_out1
	
	output:
        stdout into chunk_splitter_out

	"""
	#!/usr/bin/env python3
	
	batch_id = ("$batch_file".split('.')[0]).split('_')[2]
	with open("$batch_file", 'r') as bfile, open("$interest_ids", 'r') as inheaders:
		query_line = ''
		i = 1.0
		n = 1.0
		for line in bfile:
			query_line += (line.rstrip()) + ' '
			check = i/500
			if check == n:
				i = 0.0
				print(batch_id, '_', query_line)
				query_line = ''
			i += 1.0
		if query_line != '':
			print(batch_id, '_', query_line)
	"""
}
//chunk_splitter_out.view()


real_chuncks = chunk_splitter_out.splitText( elem: 1)
//real_chuncks.view()

process species_retriever {
	container 'python:slim-buster@sha256:9dc8667397f2d26bf648967d155bff875a523e79a964cc36d8c48576ad0df6db'

	input:
	path interest_ids from input_ids
	//path batch_file from random_batch_creator_out1
	val chunks from real_chuncks

	output:
	//stdout into species_retriever_out
	file "taxid_batch_file*.txt" into species_retriever_out

	script:
	trim_chunk = chunks.trim()
	"""
	#!/usr/bin/env python3
	import urllib.parse
	import urllib.request
	
	url = 'https://www.uniprot.org/uploadlists/'
	
	trimmed = ("${trim_chunk}").rstrip()
	one_line_query = (trimmed.split(' _ '))[1]
	identifier = (trimmed.split(' _ '))[0]
	output_filename = "taxid_batch_file" + "_" + identifier + "_" + one_line_query[0:5]+ ".txt"
	
	params = {
	'from': 'ID',
	'to': 'ACC',
	'format': 'tab',
	'columns': 'id,entry_name,organism-id,organism',
	'query': one_line_query
	}
	
	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f, open(output_filename, 'w') as outfile:
		for line in f:
			to_be_written = line.decode('utf-8')
			outfile.write(to_be_written)
	"""
}
species_retriever_out.view()
*/
