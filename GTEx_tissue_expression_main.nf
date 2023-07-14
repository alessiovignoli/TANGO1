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
        log.info 'This is the list of flags the pipeline accepts:'
        log.info '--IN_IDS              mandatory field, the file containing the pairs of IDs to be analized (do the ratio). It has to be a csv file (comma separated) '
        log.info '                      where the first field is the first ID (numerator in the ratio) a comma and second ID (denominator). It should have'
        log.info '                      a pair per line and nothing else.'
        log.info '--IN_TISSUE           mandatory flag, The sample annotation file. In version v8 this is a tab separeted file (tsv). '
        log.info '                      This file contains the relationship between the samples_IDs and the tissue of extraction.'
        log.info '                      It is though to have on each line a sample_id and the info for that sample in differrent fields, in them one is '
        log.info '                      the tissue. The file is expected to have a first line as header that is going to be excluded.'
        log.info '--IN_DATA             mandatory flag, the bulk data downloaqdable from GTEx from which the IN_IDS have to be extracted.'
        log.info '                      It is assumed that each line contains all the sample values associated to only one gene/transcript ID.'
        log.info '                      Basically every line of this file has a differerent ID per line.'
        log.info '                      It also has to have an header line (specifiable which exactly) that contains the sample-IDs in it as basically column labels/names.'
        log.info '                      From this line and the file above the relations ship between value expressions and tissues arte gathered.'
        log.info '--OUT_NAME            optional flag, default false. The output will not be renamed. If a string is given to this variable'
        log.info '                      that would be the name of the output file. On top of this since the pipeline works on each pair IDs in parallel.'
        log.info "                      there will be as many output files as there are pairs if this flag is not given. They will be named as:"
        log.info '                      IN_DATA(wiothout extension).ID1-ID2'
        log.info '                      Instead if this flag is given all the pairs files will be concatenated in a single file named OUT_NAME.'
        log.info '--GTEX_DELIMITER      optional flag, default <tab>. Since in v8 of GTEx the files are tab separated this is the default,'
        log.info '                      but it can be changed if need be.'
        log.info '--ID_POS              optional field, default 0. The column position (first = 0) that identifies where the query ID is in'
        log.info '                      the IN_TISSUE file. It is used to tell the script where to look for the IDs given.'
        log.info '--TISSUE_POS          optional field, default 6. the column in the sample annotation file containing the Tissue name.'
        log.info '                      In GTEx v8 is column 7 (6) in python notation. This field follows python notation, first column is 0.'
        log.info '--SAMPLE_POS          optional field, default 0. the column in the sample annotation file containing the sampleID.'
        log.info '                      In GTEx v8 is column 1 (0) in python notation. '
        log.info '--HEADER_LINE         optional field, default 3. The very important line in IN_DATA that contains the sample_IDs present as values of the dictionary.'
        log.info '                      From this line the values found in the line of the query IDs are correlated (through  the dictionary) to the tissue. first line = 1'
        log.info '                      this HEADER_LINE contains the sample ID present in the Dictionary computed on the IN_TISSUE file,'
        log.info '                      and that the column identified by this sample ID always contains values associated with it.'
        log.info '--OUTPUT_DIR          optional field, default params.TEST_DIR. Variable found in the nextflow.config file.'
        log.info '                      This is the directory where the output file/s are saved.'
        log.info '--STORE_DIR           optional field, default params.TEST_DIR/GTEx_data/ . Internal variable for making the pipeline more fast.'
        log.info '                      It uses the functionality of nextflow of StoreDir directive.'
        log.info '                      In brief if the file dictionary computed on IN_TISSUE is alredy existing in the STORE_DIR path the execution'
        log.info '                      of the first process is skipped. Basically if the dictionary has already been computed for this data just use that.'
        log.info '--CONTAINER           optional flag, The container used to run the piepiline, python3 is the only requirement. Better not change this.'
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



// Include section

include { GTEx_pair_id_tissue_expr } from "${params.PIPES}GTEx_tissue_expression_per_pair_ids"



workflow {
        GTEx_pair_id_tissue_expr(params.IN_IDS, params.IN_TISSUE, params.IN_DATA)
        GTEx_pair_id_tissue_expr.out.outfile.view()
        GTEx_pair_id_tissue_expr.out.stout.view()                      // for debug
}

workflow.onComplete { println 'Done' }
