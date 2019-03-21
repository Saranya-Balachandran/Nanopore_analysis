#!/usr/bin/env nextflow
// Command line input: nextflow run ../Scripts/nano_next.nf -params-file ../Scripts/config.yaml




reads = Channel
	.fromPath(params.fastq)
	.map { file -> tuple(file.baseName, file) }

genome = Channel.fromPath(params.genomeref)





/*
 * Step 1. Filter low quality data
 */

process readFilt{
	
	input:
	set datasetID, file(read_fastq) from reads
	output:
	set datasetID, file("${datasetID}_filt.fastq") into nano_filt
	

	"""
	head ${read_fastq}
	cat  ${read_fastq} | NanoFilt -q 8 -l 500 > ${datasetID}_filt.fastq
	
	"""

}


/*
 * Step 2. Preprocess the data
 */

process readCleanup{
	input:
	set datasetID, file(read_fastq) from nano_filt
  
	output:
	set datasetID, file("${datasetID}_pore.fastq") into pore_clean

	"""
	porechop -i ${read_fastq} -o ${datasetID}_pore.fastq
	
	"""

}



/*
 * Step 3. Align the fastq file using minimap2, sort and index the bam files
 */
process alignmentMap {
 	publishDir params.outputdir mode: 'copy', overwrite: false
	
	input:
	file ref_seq from genome
	
        set datasetID, file(read_fastq) from pore_clean
	
	output:
	set datasetID, file("${datasetID}_minmap.bam") into bam_file
	
	set datasetID, file("${datasetID}_sorted.bam") into bam_sort_file

	set datasetID, file("${datasetID}_sorted.bam.bai") into bam_index

	"""
	
	minimap2 -t {threads} -ax map-ont -p {params.secondary_score_ratio} -N {params.maximum_secondary} ${ref_seq} ${read_fastq} > ${datasetID}_minmap.bam
	 
	samtools sort ${datasetID}_minmap.bam > ${datasetID}_sorted.bam
	
	samtools index ${datasetID}_sorted.bam
	
	"""
}



/*
 * Step 6. Call Structural variants
 */
process variantCalling {
	publishDir params.outputdir mode: 'copy', overwrite: false
	input:
	set datasetID, file(alignmentsorted_file) from bam_sort_file
	set datasetID, file(alignmentsorted_file_bai) from bam_index

	output:
	set datasetID, file("${datasetID}.vcf") into vcf
	
	"""
		
	sniffles -s 1 -m ${alignmentsorted_file} -v ${datasetID}.vcf --genotype

	"""
}


/*
 * Step 7. Call Structural variants
 */
process vcfAnalysis {
	publishDir params.outputdir mode: 'copy', overwrite: false
	input:
	set datasetID, file(vcf_file) from vcf

	output:
	file '${datasetID}_breakpoint.csv' into breakpoint
	file '${datasetID}_deletion.csv' into deletion
	
	
	"""
		
	python {Nanopore_analysis}/Scripts/vcf_analysis.py ${vcf_file} ${datasetID}_breakpoint.csv ${datasetID}_deletion.csv

	"""
}



