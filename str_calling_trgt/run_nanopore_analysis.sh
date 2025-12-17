#! /bin/bash

### Submit this Script with: sbatch <script.sh> ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition shortterm/debug/longterm:
#SBATCH --partition=longterm
#  Use so many node:
#SBATCH --nodes=1
#  Request so many cores (hard constraint):
#SBATCH -c 4
#  Request so much of memory (hard constraint):
#SBATCH --mem=400GB
#  Find your job easier with a name:
#SBATCH --job-name=nanopore
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default


# Load your necessary modules:
module load singularity

module load nextflow/v22.04.1
#module load vep/105 # This worked before the April 2023 cluster upgrade
module load vep/105

mkdir -p $WORK/nanopore_launchdir
rsync -rc --update * $WORK/nanopore_launchdir/
cd $WORK/nanopore_launchdir
chmod u+x *
rm slurm*

# if -resume option not present, then clean the nextflow
if [[ $1 = -resume ]]; then
    echo "Resuming the previous nextflow run"
else
    echo "Cleaning up all the previous nextflow runs"
    echo "If you wanted to resume the previous run, use the '-resume' option"
    while [[ $(nextflow log | wc -l) -gt 1 ]]; do
        nextflow clean -f
    done
fi

# Submit the Nextflow Script:



nextflow run nano_next.nf -params-file config.yaml -profile omics -resume
