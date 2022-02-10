#!/usr/bin/env nextflow
// Made by: Igor Trujnara, Alessio Vignoli

nextflow.enable.dsl = 2


params.help = false


// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
	log.info "This pipeline computes kmers of length specified by the --kmer flag. "
    log.info "The input is two types of files: one FASTA file given with the --fasta option and one Phobius prediction file given with --pred."
    log.info "The prediction file is a direct stdout redirection of Phobius with -short."
    log.info "This means that for each FASTA file each sequence is found in the associated prediction file."
    log.info "The order of FASTA and prediction files does not matter."
    log.info ""
    log.info "--pairing flag, when given, no matter the content, will always use the non-canonical pairing algorithm (explained afterwards)."
    log.info "--pairing \"standard\" is the only case when the flag is explicitly mentioned, and in this case the canonical pairing algorithm is used."
    log.info "Not mentioning the pairing flag will use the canonical pairing scheme."
    log.info "Canonical pairing scheme:"
    log.info "The glob pattern in the --pred flag has to encapsulate also the FASTA files."
    log.info "e.g. predfile = sample1.txt, sample2.txt  fastafile = sample1.fa, sample2.fa"
    log.info "In this case the --pred flag will be --pred \"*.txt\" "
    log.info "Since the names are identical, the pairs generated will be sample1.txt with sample1.fa etc."
    log.info "WARNING: .txt and .fa must be used in this case. File names cannot contain other dots."
    log.info "Non-canonical pairing scheme:"
    log.info "In this case the pred file and the FASTA file name can be whatever. The only thing that links them is the content matched by the asterisk."
    log.info "e.g. predfile = sample1.txt, sample2.txt  fastafile = sample_1.fasta, sample_2.fasta"
    log.info "In this case, the flags will be --pred \"sample*.txt\" --fasta \"sample_*.fasta\""
    log.info "WARNING: Only one asterisk is recognized as glob pattern."
    log.info ""
    log.info "The last flag is --feature. It controls which type of phobius prediction feature is extracted. Possible ones are:"
    log.info "default is s for special helix, all the possible are: [c, i, o, -, n, s, l]  c = signal peptide, i = inside membrane(cytoplasm),"
    log.info "o = outside membrane, - = helix (in phobius originalmodel), (only in phobius-M7or later) => -n- = normal-helix"
    log.info "s- = special-helix and -l- = loop-inramembrane"
    log.info "they have to be given to the pipeline in one letter code"
}



// params section

params.kmer = 3
params.feature = "s"
params.pairing = "standard"
params.fasta = "data/sequence.fa"
params.pred = "data/prediction.txt"
params.outdir = "results/"
params.CONTAINER = "python:slim-buster@sha256:fe2971bedd019d952d4458afb1fe4e222ddb810150008c1dee5a068d38bb0e43"


// include section

include { ch_pairer } from "${params.PIPES}channel_pairer.nf"
include { pairer } from "${params.PIPES}input_files_pairer.nf"

/*
process findSequences {
    input:
    path pyscript
    path pred_file
    val feature_id
    
    output:
    path "seq_ranges_*.txt", emit: ranges
	stdout emit: standardout                    // debug porpouses    
    script:
    suffix = "seq_ranges_${pred_file.getFileName()}".split('\\.')[0] + ".txt"
    """
    python3 ${pyscript} ${pred_file} ${suffix} ${feature_id}
    """
}
process extractSequences {
    input:
    path pyscript
    path ranges
    path fasta
    
    output:
    //path "seqs_*.txt"
	stdout emit: standardout                    // debug porpouses
    
    script:
    suffix = "${ranges.getFileName()}".split('\\.')[0] + "__" + "${fasta.getFileName()}".split('\\.')[0] + ".txt"
    """
	python3 ${pyscript} ${ranges} ${fasta} ${suffix}
    """
}
*/

process findAndExtractSeq {
    container params.CONTAINER
    
	input:
    path findscript
    path extrscript
    tuple val(matcher), path(predFile), path(fastaFile)
    val featureID

    output:
    path "seqs_pair_*.fa", emit: sub_seqs
	stdout emit: standardout

    script:
	suffix = "seqs_${matcher}.txt"
	prefix = "${predFile}".split("${matcher}")[0].split('\\.')[0]						// used later on at the final step for the output name
    """
    python3 ${findscript} ${predFile} ${suffix} ${featureID}
    python3 ${extrscript} ${suffix} ${fastaFile} "seqs_pair_${matcher}.fa"
	echo ${prefix}
    """
}


process countKmers {
    container params.CONTAINER

    input:
    path pyscript
    path seqs
    val kmer

    output:
    path "seq_kmers_*.txt", emit: kmers_counts
    stdout emit: standardout                    // debug porpouses

    script:
	suffix = "${seqs.getFileName()}".split("seqs_pair_")[1].split('\\.')[0] + ".txt"
    """
    python3 ${pyscript} ${seqs} "seq_kmers_${suffix}" ${kmer}
    """
}

process sumKmers {
    container params.CONTAINER
    publishDir params.outdir, mode: "move", overwrite: false
    scratch true

    input:
    path pyscript
    path kmers
	val suffix

    output:
    path "total_kmers_*", emit: tot_kemers
    stdout emit: standardout

    script:
	final_name =  "total_kmers_" + "${suffix}_".replaceAll("\n", "") + "${kmers}".split("seq_kmers_")[1].replaceAll("bubba","")
    """
    for i in `echo ${kmers}`; do cat \$i >> TMP; done
    python3 ${pyscript} TMP ${final_name}
    """
}

workflow countHelixKmers {
    take:
    kmer
    featName
    prediction
    fasta

    main:
    // script finder
    findscript = params.SCRIPTS + "find_seqs.py"
    extrscript = params.SCRIPTS + "extract_seqs.py"
    countscript = params.SCRIPTS + "kmer_count.py"
    sumscript = params.SCRIPTS + "sum_kmers.py"

    // channel definitions
    pairedInputs = ""
    if(params.pairing == "standard"){
		pairedInputs = Channel.fromFilePairs( prediction + ".{txt,fa}").map{ [it[0], it[1][1], it[1][0]] }	// not to have list of lists and fa comes before txt
    } else {
		pairer(prediction, fasta)
        pairedInputs = pairer.out.right_pairs
    }
	findAndExtractSeq(findscript, extrscript, pairedInputs, featName)
    countKmers(countscript, findAndExtractSeq.out.sub_seqs, kmer)
	sumKmers(sumscript, countKmers.out.kmers_counts.collect(), findAndExtractSeq.out.standardout.last())
    
	emit:
    final_out = sumKmers.out.tot_kemers
	//stdout = sumKmers.out.standardout // for debug porpouses
}

workflow {
    countHelixKmers(params.kmer, params.feature, params.pred, params.fasta)
	countHelixKmers.out.final_out.view()
	//countHelixKmers.out.stdout.view()		// for debug porpouses
}