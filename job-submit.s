#!/bin/bash
#SBATCH --job-name=rna-seq-nf
#SBATCH --partition=p_kongt_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=94GB
#SBATCH --time=24:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
module use /projects/community/modulefiles
module load nextflow 
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"
mkdir -p $NF_Work_Dir
nextflow run rna-seq-genecount.nf --SingleEnd="true" --reads="/scratch/rw409/oarc/skin/*.fastq.gz" -w $NF_Work_Dir -resume 

