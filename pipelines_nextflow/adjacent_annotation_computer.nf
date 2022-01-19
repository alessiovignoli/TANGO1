#!/usr/bin/env nextflow

/*
*
*/

nextflow.enable.dsl=2

params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
        log.info "This pipeline computes first finds the -1 0 +1 adjacent annotations to a prediction of phobius"
        log.info "then once this is computed for each input protein the frequency of each adjecent type is computed and stored in three "
        log.info "final files with extension .domain_freq more on this files later"
        log.info ""
        log.info "the firdt input is the domain_info file that containd the uniprot information of a given protein in one line with the format:"
        log.info "	PI00068AA065 [IPR036514;sgnh hydrolase superfamily;471;634] [IPR002656;acyltransferase 3 domain;5;321] [IPR043968;sl gnh domain;400;626]"
        log.info "where the first field is the protein id and then there are as many lists as there have been found annotations for this entry in uniprot"
        log.info "there is the keyword associted with each feature a word descriptive feature and"
        log.info "and the extremities of the feature (start and end), there can be more of them"
	log.info "this file is given to the pipeline using the flag --INPUT_DOMAIN_INFO and is generated using the script:"
        log.info "	uniprot_rest_query.py    look into that for more details"
        log.info ""
        log.info "the second input file --INPUT_TXT flag is a phobius short prediction redirection file (just a > after the vommand)"
        log.info "it can be created using:	phobius_parallelization.nf	pipeline"
        log.info "but it basically stores the information about the inquired feature for example a special helix prediction"
        log.info "in thid file the ids of the domain_info are searched and when found the position of the prediction"
        log.info "are extracted (start and end of the predicted feature ex. special helix)"
        log.info "based on this info the most closest to the left (aka -1 adjacent) annnotationd are found"
        log.info "as well as the center (completely inglobing the feature predicted) aka 0 and closest to the right aka +1"
        log.info "more info in 		adjacent_annotation_finder.py"
        log.info ""
        log.info "the third and optial flag is --FIELD that is a one letter code for all the possible features predicted by phobius"
        log.info 'default is s for special helix, all the possible are: [c, i, o, -, n, s, l]  c = signal peptide, i = inside membrane(cytoplasm),'
        log.info 'o = outside membrane, - = helix (in phobius originalmodel), (only in phobius-M7or later) => -n- = normal-helix'
        log.info '-s- = special-helix and -l- = loop-inramembrane'
        log.info 'they have to be given to the pipeline in one letter code'
        log.info ""
        log.info "as mentioned the pipeline creates a file where all for example left annotation found are stored as well as the "
        log.info "start and end of both adjacent anotation and predicted feature (in this order) the location of the file is printed on screen"
        log.info "and the name is always sample_left.annot sample_center.annot sample_right.annot"
        log.info "this file comprise the summary of all that has been found and are tipically under the tmp/ dir in work/ "
        log.info "the content has the followinfg format:"
        log.info "	35 405 300 320 tpr repeat"
	log.info ""
        log.info "the name of the final three file are INPUT_DOMAIN_INFO pattern before the * or the . (if * is not present)"
        log.info "followed by _left (or center or right) and .domain_freq"
        exit 1
}

params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"
params.INPUT_DOMAIN_INFO = "${params.TEST_DIR}bubbabubba.domain_info"
params.INPUT_TXT = "${params.TEST_DIR}bubbabubba.txt"
params.OUTPUT_DIR = "${params.TEST_DIR}"
params.FIELD = "s"

// Modules dependencie section

include { pairer } from "${params.PIPES}input_files_pairer"

process  adjacency_finder {
	//publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${in_pred}" }

	input:
	tuple val(id), path(in_domain_info), path(in_pred)
	path py_script

	output:
	//stdout emit: standardout
	path "${out_name}_left.adj", emit: intermidiate_left_tmp
	path "${out_name}_center.adj", emit: intermidiate_center_tmp
	path "${out_name}_right.adj", emit: intermidiate_right_tmp

	script:
	out_name = "${in_domain_info}".split('\\.')[0]
	"""
	./${py_script} ${in_domain_info} ${in_pred} ${params.FIELD} ${out_name}
	"""
}


process summarizer_of_files {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
	container params.CONTAINER
	tag { "${all_one_type_annotations}" }

	input:
	path all_one_type_annotations
	path pyscript2

	output:
	path "${prefix}", emit: final_out
	stdout emit: standardout

	script:
	prefix = "${params.INPUT_DOMAIN_INFO}".split('/')[-1].split('\\.')[0].split('\\*')[0] + "_" +"${all_one_type_annotations}".split('_')[1].split('\\.')[0] + ".domain_freq"
	//              the above line creates the  output files based on the common part to the glob pattern that uses asterisc
	"""
	./${pyscript2} ${all_one_type_annotations} tmp
	sort -k1 -nr tmp > ${prefix}
	rm tmp
	"""
}



/*
process summarizer_of_files {
	publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
        container params.CONTAINER
        tag { "${all_one_type_annotations}" }

	input:
	path all_one_type_annotations
	
	output:
	path "${prefix}", emit: final_out
	stdout emit: standardout

	script:
	prefix = "${params.INPUT_DOMAIN_INFO}".split('/')[-1].split('\\.')[0].split('\\*')[0] + "_" +"${all_one_type_annotations}".split('_')[1].split('\\.')[0] + ".domain_freq"
	//		the above line creates the  output files based on the common part to the glob pattern that uses asterisc
	"""
	#!/usr/bin/env python3
	
	list_of_freq = []
	with open("${all_one_type_annotations}", 'r') as infile:
		for line in infile:
			domain_name = '' 	# annotation or domain string secription or name it needs the new linw character for the write
			for wordz in line.split(' ')[4:]:
				domain_name += (' ' + wordz)
			if domain_name in list_of_freq:
				i = list_of_freq.index(domain_name)
				list_of_freq[(i+1)] += 1
			else:
				list_of_freq.append(domain_name)
				list_of_freq.append(1)
	with open("${prefix}", 'w') as freqfile:
		for n, annotation in enumerate(list_of_freq[::2]):		# proceding two by two   n will not represent the real index of annotation
			freqfile.write(( str(list_of_freq[(n*2+1)]) + annotation))
			#print(str(list_of_freq[(n*2+1)]) + annotation, end='')
	"""
}
*/


workflow adjacent_domains_compiler {

	take:
	pattern_to_domaininfo
	patter_to_shortpred

	main:
	pairer(pattern_to_domaininfo, patter_to_shortpred)
	adj_finder_pyscript = params.SCRIPTS + "adjacent_annotation_finder.py"
	adjacency_finder(pairer.out.right_pairs, adj_finder_pyscript)
	adjacency_finder.out.intermidiate_left_tmp.collectFile(name: "sample_left.annot" ).set{ tmp_left }
	adjacency_finder.out.intermidiate_center_tmp.collectFile(name: "sample_center.annot" ).set{ tmp_center }
	adjacency_finder.out.intermidiate_right_tmp.collectFile(name: "sample_right.annot" ).set{ tmp_right }
	tmp_left.concat( tmp_center, tmp_right ).set{ proper_tmp }
	adj_freq_compute_pyscript = params.SCRIPTS + "uniprot_annotation_domain_processer.py"
	summarizer_of_files(proper_tmp, adj_freq_compute_pyscript)

	emit:
	intermidiate_exhaustive = proper_tmp
	stout = summarizer_of_files.out.final_out
}

workflow {
	adjacent_domains_compiler(params.INPUT_DOMAIN_INFO, params.INPUT_TXT)
	adjacent_domains_compiler.out.intermidiate_exhaustive.view()
	adjacent_domains_compiler.out.stout.view()
}

