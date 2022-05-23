#!/usr/bin/env nextflow

OMA_folder = "/home/alessio/Desktop/erasmus-internship/OMA/"
in_fasta = Channel.fromPath(OMA_folder + "*_raw.fasta")

process phobius_image_test {
	container 'phobius_image'

	input:
	path fastas from in_fasta

	output:
	stdout phobius_image_test_out

	script:
	"""
	perl /home/phobius101_linux/tmp/tmpbyom7j/phobius/phobius.pl ${fastas}
	"""
}
phobius_image_test_out.view()
