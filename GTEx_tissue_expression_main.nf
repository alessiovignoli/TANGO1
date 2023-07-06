#!/usr/bin/env nextflow


// If the --help parameter is used in the command line, the pipeline will print the following help section:

if (params.help) {
        log.info 'This is the help section of this pipeline'
        log.info 'This pipeline wiull compute an average count per tissue of the requested transcripts/gene/exon ecc.. ids on GTEx data.'
        log.info 'It has been written for GTEx v8. But in theory all that is needed is a valid ID corresponding to the respective data file.'
        log.info 'In GTex data there are various type of publicly available RNA-seq data: gene count, transcripts counts, exon counts ecc..'
        log.info 'The structure of this type of file as this version is the following:'
        log.info 'Header line with sample ids'
        log.info 'id1_col maybe_id1_col2 sample1_val sample2_val sample3_val .... sample17234_val'
        log.info 'for many lines usually gzipped.'
        log.info ''
        log.info 'Since the number of samples is so big it is hard to make sense of it. This pipeline gathers the values in the sample columns'
        log.info 'based on the tissue type of that column sample_id. Since the tissue types are in the order of 60 it becames more understandable.'
        log.info 'This is done per transcript/gene or genarally speacking ID. So that at the end it will look something like:'
        log.info 'requested_ID Brain-frontal_average Muscle-skeletal_average ecc..'
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        log.info ''
        exit 1
}


params.CONTAINER = "python@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" // Python 3.9.5 slim buster
params.IN_IDS = null
params.IN_TISSUE = null
params.IN_DATA = null
params.OUT_NAME = false
params.GTEX_DELIMITER = '\t'
params.ID_POS = 0
params.TISSUE_POS = 6
params.SAMPLE_POS = 0
params.HEADER_LINE = 3
params.OUTPUT_DIR = "${params.TEST_DIR}"	
params.STORE_DIR = "${params.TEST_DIR}GTEx_data/"
params.PUBLISH = true
params.SCRATCH = true


// Include section

include { GTEx_id_tissue_expr } from "${params.PIPES}GTEx_tissue_expression_per_id"



workflow {
        GTEx_id_tissue_expr(params.IN_IDS, params.IN_TISSUE, params.IN_DATA)
        GTEx_id_tissue_expr.out.outfile.view()
        GTEx_id_tissue_expr.out.stout.view()                      // for debug
}

workflow.onComplete { println 'Done' }
