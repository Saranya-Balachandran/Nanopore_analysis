#!/usr/bin/env nextflow
// Command line input: nextflow run /pathtoscript/nano_next.nf --fastq 'merged_fastq_files' --ref 'reference.fa'


params.fastq = Channel.fromPath( 'fastq_runid.fastq')
params.ref = Channel.fromPath( 'fasta_runid.fasta')
reads = file(params.fastq)
genome = file(params.ref)

/*
 * Step 1. Filter low quality data
 */

process readFilt{
	input:
	file read_fastq from reads
  
	output:
	file 'filt.fastq' into nano_filt

	"""
	cat $read_fastq | NanoFilt -q 8 -l 500 > filt.fastq
	
	"""

}


/*
 * Step 2. Preprocess the data
 */

process readCleanup{
	input:
	file read_fastq from nano_filt
  
	output:
	file 'read.fastq' into pore_clean

	"""
	porechop -i ${read_fastq} -o read.fastq
	
	"""

}



/*
 * Step 3. Align the fastq file using minimap2, sort and index the bam files
 */
process alignmentMap {

	input:
	file ref_seq from genome
	
        file read_fastq from pore_clean
	
	output:
	file 'read.bam' into bam_file

	"""
	
	ngmlr -r  ${ref_seq} -q ${read_fastq} -o read.bam

	
	"""
}

/*
 * Step 4. sort the bam files
 */
process alignmentSort {

	input:
	file alignment_file from bam_file
	
	
	output:
	file 'reads_sorted.bam' into bam_sort_file

	"""
		
	samtools sort ${alignment_file} > reads_sorted.bam
	
	"""
}
/*
 * Step 5. index the bam files
 */
process alignmentIndex {

	input:
	file alignmentsorted_file from bam_sort_file
	
	
	
	
	"""
		
	samtools index ${alignmentsorted_file}

	"""
}

/*
 * Step 6. Call Structural variants
 */
process variantCalling {

	input:
	file alignmentsorted_file from bam_sort_file

	
	
	
	"""
		
	sniffles -m ${alignmentsorted_file} -v final.vcf

	"""
}






result.subscribe { println it }

