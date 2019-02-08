#!/usr/bin/env nextflow
// Command line input: nextflow run /pathtoscript/stats_check.nf --in 'merged_fastq_files'


params.in = Channel.fromPath( 'fastq_runid.fastq')
seq = file(params.in)


process runStats {

	input:
	file 'input.fastq' from seq
	
	

	"""
		
	NanoStat --fastq input.fastq -n Statsreport
	
	NanoPlot -t 2 --fastq input.fastq --maxlength 40000 --plots hex dot

	"""
}

result.subscribe { println it }

