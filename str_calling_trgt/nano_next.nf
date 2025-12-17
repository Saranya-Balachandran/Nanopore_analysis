#!/usr/bin/env nextflow
// Command line input: nextflow run /home/nanopore/Saranya/Scripts/nano_next.nf -params-path(/home/nanopore/Saranya/Scripts/config.yaml

//plot_str_trgt(call_str_trgt.out.str_vcf,call_str_trgt.out.repeat_bed,genome, annotate_repeat_expansions.out.repeat_ann)
nextflow.enable.dsl=2

/*Global parameters*/
reads = channel
	.fromPath(params.fastq)
	.map {it -> [it.baseName, it]}.groupTuple()

genome = channel.fromPath(params.genomeref,checkIfExists: true).first()

str_list = "${launchDir}/data/wf_str_repeats.bed"
variant_catalogue_hg38 = "${launchDir}/data/variant_catalog_hg38.json"

trid = Channel
    .fromPath("${launchDir}/data/wf_str_repeats.bed")
    .splitCsv(sep: '\t')
    .map { row -> (row[4]) }


/*---------------Main workflow-------------*/
workflow{
readFilt(reads)

readCleanup(readFilt.out.nano_filt)
genomeindex(genome)
alignment_minimap(readCleanup.out.pore_clean, genomeindex.out.index, genome)
call_str_trgt(alignment_minimap.out.minimap_bam_sort_path,alignment_minimap.out.minimap_bam_index, genome, str_list)
annotate_repeat_expansions(call_str_trgt.out.str_vcf, variant_catalogue_hg38)
call_str_trgt.out.str_vcf
   .combine(trid)
   .set {channel_plot}
plot_str_trgt(channel_plot,genome)
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
	/opt/conda/envs/nanopore/bin/minimap2 -t ${params.threads} -K 500M -r 50000,50000 -end-bonus=10000 -no-end-flt -ax map-ont ${ref_seq} ${read_fastq} > ${datasetID}_minimap.bam

	/opt/conda/envs/nanopore/bin/samtools sort -@ ${params.threads} ${datasetID}_minimap.bam > ${datasetID}_minimap_sorted.bam

	/opt/conda/envs/nanopore/bin/samtools calmd -@ ${params.threads} -b ${datasetID}_minimap_sorted.bam ${ref} > ${datasetID}_minimap_sorted_calmd.bam

	/opt/conda/envs/nanopore/bin/samtools index -@ ${params.threads} ${datasetID}_minimap_sorted_calmd.bam

	"""
}

/*
* STEP5 TRGT_str caller
*/
process call_str_trgt {

	input:
  	tuple val(datasetID), path(bam)
	tuple val(datasetID), path(bam_index)
	path(ref)
	path(repeat_bed)
	output:
	tuple val(datasetID), path("*.vcf.gz*"), path("*.spanning.bam*"),path("repeats_subset.bed"), emit: str_vcf
	"""
	awk -F '\\t' 'BEGIN{OFS=FS} {m=\$4; i=\$5; \$4="MOTIFS=" m ";ID=" i ";STRUC=<TR>"; print}' ${repeat_bed} | awk -F '\\t' 'NF{NF-=2};1'| tr -d \$'\r' > repeats_subset.bed
	if [[ -s repeats_subset.bed ]]; then
	samtools faidx ${ref}
	head repeats_subset.bed
	trgt genotype --genome ${ref} \
	--repeats repeats_subset.bed \
	--reads ${bam} \
	--max-depth 50000 \
	--output-prefix ${datasetID}_str

	bcftools sort ${datasetID}_str.vcf.gz -o ${datasetID}_str.vcf.gz -O b

	tabix -p vcf ${datasetID}_str.vcf.gz
		else
	            echo "blank subset BED"
	        fi


	"""


}

process plot_str_trgt{
	errorStrategy = 'ignore'
	input:
			tuple val(datasetID), path(vcf), path(bam),path(repeat_bed), val(trid)
			path(ref)
	output:
			tuple val(datasetID), path("*.spanning.sorted.bam*"), path("*.svg")

	"""
		samtools sort -o ${datasetID}_str.spanning.sorted.bam ${datasetID}_str.spanning.bam
		samtools index ${datasetID}_str.spanning.sorted.bam

	samtools faidx ${ref}
	trgt plot --genome ${ref} \
		  --repeats repeats_subset.bed \
		  --vcf ${datasetID}_str.vcf.gz \
		  --spanning-reads ${datasetID}_str.spanning.sorted.bam \
		  --repeat-id ${trid} \
		  --max-allele-reads 15 \
		  --image ${datasetID}.${trid}.svg

	"""
}
process annotate_repeat_expansions {
    // annotate using Stranger
    input:
        tuple val(datasetID), path(vcf), path(bam), path(repeats_bed)
        path(variant_catalogue_hg38)
    output:
        tuple val(datasetID), path("${datasetID}_repeat-expansion_annotated.vcf.gz"),path("${datasetID}_repeat-expansion_annotated.vcf.gz.tbi")
	path("*_annotated_filtered.tsv"), emit: repeat_ann
    script:
        """
        stranger -f ${variant_catalogue_hg38} -t ${datasetID}_str.vcf.gz \
            | sed 's/\\ /_/g' \
            | bgzip -c > ${datasetID}_repeat-expansion_annotated.vcf.gz
        tabix -p vcf ${datasetID}_repeat-expansion_annotated.vcf.gz
        java -jar /opt/snpEff/SnpSift.jar extractFields ${datasetID}_repeat-expansion_annotated.vcf.gz \
            CHROM POS ALT REF GEN[*].SD GEN[*].MC MOTIFS TRID STR_NORMAL_MAX STR_PATHOLOGIC_MIN Disease STR_STATUS > ${datasetID}_repeat-expansion_annotated.tsv
        sed -i -e "1s|GEN\\[\\*\\]\\.SD|#spanning_read|" -e "1s|GEN\\[\\*\\]\\.MC|#repeat_length|" ${datasetID}_repeat-expansion_annotated.tsv
	awk 'FNR==1 {print \$NF,\$0}' ${datasetID}_repeat-expansion_annotated.tsv > ${datasetID}_repeat-expansion_annotated_filtered.tsv	
	awk 'NR==1 || FNR>1 {if(\$5 > 20){print}}' ${datasetID}_repeat-expansion_annotated.tsv >> ${datasetID}_repeat-expansion_annotated_filtered.tsv
	"""
}
