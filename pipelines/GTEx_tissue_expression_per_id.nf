#!/usr/bin/env nextflow



process group_sample_per_tissue {
	container params.CONTAINER

	input:
	path tissue_info

	output:
	stdout emit: standardout                                              // for debug

	script:
        """
	gtex_group_sample_per_tissue.py --sample_annotations ${tissue_info} \
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
	stout = group_sample_per_tissue.out.standardout			//for debug

}


