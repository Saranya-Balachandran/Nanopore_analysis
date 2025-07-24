#!/bin/bash

#SBATCH --mem=128GB
#SBATCH --partition=shortterm
#SBATCH --cpus-per-task=32

## https://nanoporetech.github.io/modkit/intro_pileup.html #description-of-bedmethyl-output
module load singularity
module load nextflow
export NXF_SINGULARITY_HOME_MOUNT=true

cd $SCRATCH

nextflow run /data/humangen_external/Nanopore/wf-human-variation \
	--bam '/data/humangen_external/Daniela/MakuP5Monarch/bam_pass/' \
	--ref '/data/humangen_external/Nanopore/reference/hg38_ebv.fa' \
	--sample_name 'MakuP5Monarch' \
	--out_dir "/data/humangen_external/Daniela/MakuP5Monarch/"
	--sv \
	--mod \
	--phased \
	--bam_min_coverage 8 \
	--override_basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.3.0' \
	-profile singularity,uzl_omics

