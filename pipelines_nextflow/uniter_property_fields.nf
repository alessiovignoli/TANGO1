#!/usr/bin/env nextflow

/*
*  This module is created for solving the issue with re use of the process properties_combiner
*  it is not written to be executed by itself, this module is part of the phisic_properties_combiner.nf workflow
*/

nextflow.enable.dsl=2


params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43" 

process properties_combiner {
        //publishDir(params.OUTPUT_DIR, mode: 'copy', overwrite: false)
        tag { "${output_name}" }
        container params.CONTAINER

        input:
        path hydro_file
        path average_file
        path aacomp_file
        path pyscript1

        output:
        path "${output_name}", emit: outname
        stdout emit: standardout

        script:
	output_name = "${hydro_file}".split('\\.')[0] + ".${params.SEC_KEYWORD}_uni"
        """
        ./${pyscript1} ${hydro_file} ${average_file} ${aacomp_file} 2>tmp.err 1>${output_name}
	cat tmp.err
        """
}

workflow phisic_prop_louncher {

        take:
        hydro
        average
        aacomp

        main:
        phisic_properties_py = params.SCRIPTS + "properties_combiner.py"
        properties_combiner(hydro, average, aacomp, phisic_properties_py)

	emit:
        stout1 = properties_combiner.out.standardout
	unified_phi =  properties_combiner.out.outname
	
}

