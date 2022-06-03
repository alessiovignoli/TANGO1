




laDir="${baseDir}"
//params.output="${leilaDir}/"
params.output = "${params.TEST_DIR}"
params.aln="${leilaDir}/aln_fa/PF*"
params.gammaRate="1.0"
params.seedValue="5"
params.replicatesNum="100"
params.trimmer="trimal"
params.mode="8"
params.CONTAINER = 'lmansouri/phylo_3d:220225'

if ( params.aln ) {
  Channel
  .fromPath(params.aln)
  .map { item -> [ item.baseName.replace("_tmalign","") , item] }
  .into { fasta_aln ; original_msa_ch; oriSeqs}
}


process blocking_aln {
  tag"${id}"
  publishDir "${params.output}/trimal_automated1_trimmed_alns", mode: 'copy', overwrite: true
  input:
  set val(id), file(fasta) from fasta_aln
  
  output:
  set val(id), file("*.ph") into trimmed_aln, blocked_aln
  set val(id), file("*.fa") into trimmed_fasta, blocked_fasta
  script:
  """
  trimal -in ${fasta} -out ${id}_trimmal_aln.ph -phylip -automated1
  trimal -in ${fasta} -out ${id}_trimmal_aln.fa -automated1
  """
}

process computing_ME_trees{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/FASTME_trees", mode: 'copy', overwrite: true
  
  input:
  set val(id), file(phylip) from trimmed_aln

  output:
  set val(id), file("*.nwk"), file("*.replicates") into fastme_ch

  script:
  """
  fastme -i ${phylip} -o ${id}_tmalign_FASTME.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b ${params.replicatesNum} -B ${id}_tmalign.FastMEtree.replicates
  """
}

process computing_ML_trees{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/ML_trees", mode: 'copy', overwrite: true
  
  input:
  set val(id), file(phylip) from blocked_aln

  output:
  set val(id), file("*.treefile"), file("*.boottrees") into iqtree_ch

  script:
  """
  iqtree -s ${phylip} -b ${params.replicatesNum} 
  """
}


workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
