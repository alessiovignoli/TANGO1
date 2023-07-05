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
					--delimiter ${params.TISSUE_DELIMITER}
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

	emit:
	outfile = group_sample_per_tissue.out.tissue_sample_dict
	stout = 'bubba' //group_sample_per_tissue.out.standardout			//for debug

}


