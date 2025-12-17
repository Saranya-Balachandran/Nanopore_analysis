#!/bin/bash

#SBATCH --mem=128GB
#SBATCH --partition=shortterm
#SBATCH --cpus-per-task=32

## https://nanoporetech.github.io/modkit/intro_pileup.html #description-of-bedmethyl-output
module load singularity
module load nextflow
export NXF_SINGULARITY_HOME_MOUNT=true
export SINGULARITY_CACHEDIR="$WORK/singularity"

singularity pull docker://saranyabalachandr/dorado_1.0.4_hac_v5.2:v1

singularity exec --nv docker://saranyabalachandr/dorado_1.0.4_hac_v5.2:v1 /dorado-1.0.4-linux-x64/bin/dorado basecaller "/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v5.2.0" . --device auto --min-qscore 15 --recursive --reference /data/humangen_external/Nanopore/reference/hg38_ebv.fa --max-reads 100000000 --kit-name "SQK-LSK114" --trim adapters /Volumes/Humangenetik/Data/Nanopore/Daniela/BL-CL/MakuP5Monarch/20241113_1056_P2S-01141-B_PBA60521_7fe828ca/pod5/ > MakuP5Monarch.bam
# for SUP models run below
#singularity exec --nvdocker pull saranyabalachandr/dorado_0.8_sup:v1 /dorado-0.8.1-linux-x64/bin/dorado basecaller "/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0" . --device auto --min-qscore 15 --recursive --reference /data/humangen_external/Nanopore/reference/hg38_ebv.fa --max-reads 100000000 --kit-name "SQK-LSK114" --trim adapters /Volumes/Humangenetik/Data/Nanopore/Daniela/BL-CL/MakuP5Monarch/20241113_1056_P2S-01141-B_PBA60521_7fe828ca/pod5/ > MakuP5Monarch.bam
