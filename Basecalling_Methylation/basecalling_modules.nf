/*
 * Step 1. Basecalling
 */

process basecalling {
	label "gpu"
	publishDir "${outfolder}", mode: 'copy', overwrite: true
	input:
	tuple val(unique_name), path(pod5_path), val(outfolder)
	output:
	tuple val(unique_name), path("*.bam"), val(outfolder)

	"""
	/dorado-1.4.0-linux-x64/bin/dorado basecaller ${params.dorodo_model} --device auto --min-qscore 15 --recursive --reference ${params.reference} --max-reads ${params.max_reads} --kit-name ${params.kit_name} --trim adapters ${pod5_path} > ${unique_name}.bam
       """

}

/*
 * Step 2. Dumux to fastq
 * dorado-1.4.0-linux-x64/bin/dorado demux --threads 16 --sort-bam --no-trim --kit-name ${params.kit_name} --output-dir ${outfolder} --emit-fastq ${bam}
 */

process demux {
	label "cpu"
	publishDir "${outfolder}", mode: 'copy', overwrite: true
	input:
	tuple val(unique_name), path(bam), val(outfolder)
	output:
	tuple val(unique_name), path("*.fastq")

	"""
        samtools fastq  -T MM,ML ${bam} > ${unique_name}.fastq
        rsync -av ${unique_name}.fastq ${outfolder}	
	"""

}
