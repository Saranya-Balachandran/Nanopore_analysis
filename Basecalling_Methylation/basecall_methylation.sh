#creation_date:06.03.2026

#!/bin/bash

#SBATCH --mem=490000
#SBATCH --partition=longterm
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1

## https://nanoporetech.github.io/modkit/intro_pileup.html #description-of-bedmethyl-output
module load singularity
module load nextflow
module load nvidia-cuda
export NXF_SINGULARITY_HOME_MOUNT=true
export SINGULARITY_CACHEDIR="$WORK/singularity"


#replace the model based on https://software-docs.nanoporetech.com/dorado/latest/models/list/

singularity exec --nv docker://saranyabalachandr/dorado_1.4.0:v1 /dorado-1.4.0-linux-x64/bin/dorado basecaller "/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v5.2.0" . --device auto --min-qscore 15 --recursive --reference /data/humangen_external/Nanopore/reference/hg38_ebv.fa --max-reads 100000000 --kit-name "SQK-LSK114" --trim adapters /Volumes/Humangenetik/Data/Nanopore/Daniela/BL-CL/MakuP5Monarch/20241113_1056_P2S-01141-B_PBA60521_7fe828ca/pod5/ > MakuP5Monarch.bam
