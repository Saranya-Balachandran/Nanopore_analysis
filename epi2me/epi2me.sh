#!/bin/bash

#SBATCH --mem=128GB
#SBATCH --partition=shortterm
#SBATCH --cpus-per-task=32

module load singularity
module load nextflow
export NXF_SINGULARITY_HOME_MOUNT=true
export SINGULARITY_CACHEDIR="$WORK/singularity"

#cp -r wf-human-variation $WORK/Nanopore_workdir/
cd $WORK/Nanopore_workdir
bam_file="$1"
sample_name="$2"
out_dir="$3"

# Documenting the input parameters
echo "-------------------------------------"
echo "Message: These are the run parameters"
echo "Sample: ${sample_name}"
echo "Bam files: ${bam_file}"
echo "Out dir: ${out_dir}"
echo "-------------------------------------"

nextflow -Duser.country=US -Duser.language=en run epi2me-labs/wf-human-variation \
	--bam $bam_file \
	--ref '/data/humangen_external/Nanopore/reference/hg38.fa' \
	--sample_name $sample_name \
	--snp \
	--sv \
	--mod \
	--str \
	--phased \
	--override_basecaller_cfg 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0' \
	--bam_min_coverage 0 \
	--out_dir $out_dir \
	-profile singularity
