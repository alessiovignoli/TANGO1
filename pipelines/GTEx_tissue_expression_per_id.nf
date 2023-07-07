#!/usr/bin/env nextflow



process group_sample_per_tissue {
	container params.CONTAINER
	storeDir params.STORE_DIR

	input:
	path tissue_info

	output:
	path tissue_dict_out, emit: tissue_sample_dict
	//stdout emit: standardout                                              // for debug

	script:
	tissue_dict_out = "${tissue_info.baseName}" + ".pkl"
        """
	gtex_group_sample_per_tissue.py --sample_annotations ${tissue_info} \
					--out_name ${tissue_dict_out} \
					--tissue_pos ${params.TISSUE_POS} \
					--sample_pos ${params.SAMPLE_POS} \
					--delimiter ${params.GTEX_DELIMITER}
	"""
}


process id_tissue_average {
	container params.CONTAINER

	input:
        path tissue_info
	path bulk_data
	each ID

	output:
	path out_name, emit: id_tissue_avg_expr
	stdout emit: standardout                                              // for debug

	script:
	out_name = "${bulk_data.baseName}" + ".${ID}"
	"""
	gtex_id_expr_per_tissue.py --tissue_dict ${tissue_info} \
				   --gtex_data  ${bulk_data} \
				   --ID ${ID} \
 				   --out_name ${out_name} \
				   --id_pos ${params.ID_POS} \
				   --delimiter ${params.GTEX_DELIMITER} \
				   --header_line ${params.HEADER_LINE}
	"""
}



workflow GTEx_id_tissue_expr {

	take:
	path_to_ids
	path_to_tissue_info
	path_to_bulk_data

	main:
	in_ids =  Channel.fromPath(path_to_ids)
	in_tissue_info = Channel.fromPath(path_to_tissue_info)
	in_bulk_data = Channel.fromPath(path_to_bulk_data)

	// Build the dictionary that connects sample ids to tissue
	group_sample_per_tissue(in_tissue_info)

	// Extract IDs to work one ID at time
	in_ids.splitText(){ it.trim() }.set{ all_IDS }
	id_tissue_average(group_sample_per_tissue.out.tissue_sample_dict, in_bulk_data, all_IDS)


	// if params.OUT_NAME is not false collect the output files IDs in one single file
	outfile = null
	if (params.OUT_NAME) {
		outfile = id_tissue_average.out.id_tissue_avg_expr.collectFile( name: params.OUT_NAME, storeDir: params.OUTPUT_DIR ) 
	} else {
		outfile = id_tissue_average.out.id_tissue_avg_expr.collectFile( storeDir: params.OUTPUT_DIR)
	}


	emit:
	outfile // id_tissue_average.out.id_tissue_avg_expr
	stout = id_tissue_average.out.standardout			//for debug

}


