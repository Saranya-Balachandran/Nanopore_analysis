#!/usr/bin/env nextflow
// Command line input: nextflow run /home/nanopore/Saranya/Scripts/nano_next.nf -params-path(/home/nanopore/Saranya/Scripts/config.yaml


nextflow.enable.dsl=2

/*Global parameters*/
reads = channel
	.fromPath(params.fastq)
	.map {it -> [it.baseName, it]}.groupTuple()

genome = channel.fromPath(params.genomeref,checkIfExists: true).first()




/*---------------Main workflow-------------*/
workflow{
readFilt(reads)

readCleanup(readFilt.out.nano_filt)
genomeindex(genome)
alignment_minimap(readCleanup.out.pore_clean, genomeindex.out.index, genome)
}

/*------------------------------------*/

/*
 * Step 1. Filter low quality data
 */

process readFilt{
	label "alignment"
	input:
	tuple val(datasetID), path(read_fastq)
	output:
	tuple val(datasetID), path("${datasetID}_filt.fastq.gz"), emit:nano_filt


	"""

	gunzip -c ${read_fastq} | /opt/conda/envs/nanopore/bin/NanoFilt -q 8 | gzip > ${datasetID}_filt.fastq.gz

	"""

}

/*
 * Step 2. Preprocess the data
 */

process readCleanup{
	label "alignment"
	input:
	tuple val(datasetID), path(read_fastq)
	output:
	tuple val(datasetID), path("${datasetID}_pore.fastq.gz"), emit:pore_clean

	"""
	/opt/conda/envs/nanopore/bin/porechop -i ${read_fastq} -o ${datasetID}_pore.fastq.gz
	"""

}


/*
 * Step 3. Preprocess the data
 */

process genomeindex{
	label "alignment"
	input:
	path(ref_seq)

	output:
	path('genome.mmi'), emit:index

	"""
	/opt/conda/envs/nanopore/bin/minimap2 -t {params.threads} -I 1000G -d genome.mmi ${ref_seq}

	"""

}




/*
 * Step 4. Align the fastq path(using minimap2, sort and index the bam paths
 */
process alignment_minimap {
	label "alignment"

	input:
        tuple val(datasetID), path(read_fastq)
	path(ref_seq)
	path(ref)
	output:

	tuple val(datasetID), path("${datasetID}_minimap_sorted_calmd.bam"), emit:minimap_bam_sort_path

	tuple val(datasetID), path("${datasetID}_minimap_sorted_calmd.bam.bai"), emit:minimap_bam_index



	"""
	/opt/conda/envs/nanopore/bin/minimap2 -t ${params.threads}  -ax map-ont ${ref_seq} ${read_fastq} > ${datasetID}_minimap.bam

	/opt/conda/envs/nanopore/bin/samtools sort -@ ${params.threads} ${datasetID}_minimap.bam > ${datasetID}_minimap_sorted.bam

	/opt/conda/envs/nanopore/bin/samtools calmd -@ ${params.threads} -b ${datasetID}_minimap_sorted.bam ${ref} > ${datasetID}_minimap_sorted_calmd.bam

	/opt/conda/envs/nanopore/bin/samtools index -@ ${params.threads} ${datasetID}_minimap_sorted_calmd.bam

	"""
}
