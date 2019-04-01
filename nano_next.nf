#!/usr/bin/env nextflow
// Command line input: nextflow run ../Scripts/nano_next.nf -params-file ../Scripts/config.yaml




reads = Channel
	.fromPath(params.fastq)
	.map { file -> tuple(file.baseName, file) }

genome = file(params.genomeref)


/*
 * Step 1. Filter low quality data
 */

process readFilt{
	
	input:
	set datasetID, file(read_fastq) from reads
	output:
	set datasetID, file("${datasetID}_filt.fastq") into nano_filt
	

	"""
	
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
	stdout result
	set datasetID, file("${datasetID}_pore.fastq") into pore_clean, pore_clean1

	"""
	echo ${datasetID}
	
	porechop -i ${read_fastq} -o ${datasetID}_pore.fastq
	
	"""

}


/*
 * Step 3. Preprocess the data
 */

process genomeindex{
	input:
	file ref_seq from genome
  
	output:
	file 'genome.mmi' into index

	"""
	minimap2 -t {params.threads} -I 1000G -d genome.mmi ${ref_seq}
	
	"""

}




process alignment_ngmlr {
 	publishDir params.outputdir, mode: 'copy', overwrite: false
	
	input:
		
        set datasetID, file(read_fastq) from pore_clean
	
	file ref_seq_ngmlr from genome
	output:

	set datasetID, file("${datasetID}_ngmlr.bam") into ngmlr_bam_file
	
	set datasetID, file("${datasetID}_ngmlr_sorted.bam") into ngmlr_bam_sort_file

	set datasetID, file("${datasetID}_ngmlr_sorted.bam.bai") into ngmlr_bam_index

	"""
	
	
	ngmlr -t ${params.threads} -r ${ref_seq_ngmlr} -q ${read_fastq} --skip-write -o ${datasetID}_ngmlr.bam -x ont
	 
	samtools sort -@ ${params.threads} ${datasetID}_ngmlr.bam > ${datasetID}_ngmlr_sorted.bam
	
	samtools index -@ ${params.threads} ${datasetID}_ngmlr_sorted.bam
	
	"""
}

/*
 * Step 4. Align the fastq file using minimap2, sort and index the bam files 
 */
process alignment_minimap {
 	publishDir params.outputdir, mode: 'copy', overwrite: false
	
	input:
	
	
        set datasetID, file(read_fastq) from pore_clean1
	file ref_seq from index
	file ref from genome
	output:
	set datasetID, file("${datasetID}_minmap.bam") into minimap_bam_file

	set datasetID, file("${datasetID}_minimap_sorted.bam") into minimap_bam_sort_file_initial
	
	set datasetID, file("${datasetID}_minimap_sorted_calmd.bam") into minimap_bam_sort_file

	set datasetID, file("${datasetID}_minimap_sorted_calmd.bam.bai") into minimap_bam_index
	
	

	"""
	
	minimap2 -t {params.threads} -ax map-ont -p {params.secondary_score_ratio} -N {params.maximum_secondary} ${ref_seq} ${read_fastq} > ${datasetID}_minmap.bam
		 
	samtools sort -@ ${params.threads} ${datasetID}_minmap.bam > ${datasetID}_minimap_sorted.bam
	
	samtools calmd -@ ${params.threads} -b ${datasetID}_minimap_sorted.bam ${ref} > ${datasetID}_minimap_sorted_calmd.bam

	samtools index -@ ${params.threads} ${datasetID}_minimap_sorted_calmd.bam
	
	"""
}



/*
 * Step 6. Call Structural variants
 */
process variantCalling {
	publishDir params.outputdir, mode: 'copy', overwrite: false
	input:
	set datasetID, file(alignmentsorted_file_ngmlr) from ngmlr_bam_sort_file
	set datasetID, file(alignmentsorted_file_bai_ngmlr) from ngmlr_bam_index
	set datasetID, file(alignmentsorted_file_minimap) from minimap_bam_sort_file
	set datasetID, file(alignmentsorted_file_bai_minimap) from minimap_bam_index
	
	output:
	set datasetID, file("${datasetID}_ngmlr.vcf") into vcf_ngmlr
	set datasetID, file("${datasetID}_minimap.vcf") into vcf_minimap
	
	"""
		
	sniffles -s 1 -f 0.1 --skip_parameter_estimation -m ${alignmentsorted_file_ngmlr} -v ${datasetID}_ngmlr.vcf --genotype
	sniffles -s 1 -f 0.1 --skip_parameter_estimation -m ${alignmentsorted_file_minimap} -v ${datasetID}_minimap.vcf --genotype

	"""
}


/*
 * Step 7. Call Structural variants
 */
process vcfAnalysis {
	publishDir params.outputdir, mode: 'copy', overwrite: false
	input:
	set datasetID, file(vcf_file_ngmlr) from vcf_ngmlr
	set datasetID, file(vcf_file_minimap) from vcf_minimap
	output:
	file "${datasetID}_minimap_breakpoint.csv" into vcf_break_minimap
	file "${datasetID}_minimap_deletion.csv" into vcf_del_minimap
	file "${datasetID}_ngmlr_breakpoint.csv" into vcf_break_ngmlr
	file "${datasetID}_ngmlr_deletion.csv" into vcf_del_ngmlr
	
	
	"""
		
	python /home/nanopore/Saranya/Scripts/vcf_analysis.py ${vcf_file_minimap} ${datasetID}_minimap_breakpoint.csv ${datasetID}_minimap_deletion.csv

	python /home/nanopore/Saranya/Scripts/vcf_analysis.py ${vcf_file_ngmlr} ${datasetID}_ngmlr_breakpoint.csv ${datasetID}_ngmlr_deletion.csv

	"""
}


