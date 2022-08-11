#!/bin/bash
#SBATCH --job-name=SRA_idenogs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=idengos@albany.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=sra.log

fastq-dump --split-files SRR8514718
fastq-dump --split-files SRR8514719
fastq-dump --split-files SRR8514720
fastq-dump --split-files SRR8514721
fastq-dump --split-files SRR8514722
fastq-dump --split-files SRR8514723
fastq-dump --split-files SRR8514724
fastq-dump --split-files SRR8514725
fastq-dump --split-files SRR8514726
fastq-dump --split-files SRR8514727
fastq-dump --split-files SRR8514728
fastq-dump --split-files SRR8514729
fastq-dump --split-files SRR8514730
fastq-dump --split-files SRR8514731
fastq-dump --split-files SRR8514732
fastq-dump --split-files SRR8514733
fastq-dump --split-files SRR8514734
fastq-dump --split-files SRR8514735
fastq-dump --split-files SRR8514736
fastq-dump --split-files SRR8514737
fastq-dump --split-files SRR8514738
