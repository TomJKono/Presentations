#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -A mcgaughs
#PBS -m abe
#PBS -M user@umn.edu
#PBS -q batch

# Load modules
module load bwa/0.7.15

# Find all R1 samples
SAMPLES_F=($(find /panfs/roc/scratch/konox006/reads -type f -name '*R1.fastq.gz' | sort -V))
# Find all R2 samples
SAMPLES_R=($(find /panfs/roc/scratch/konox006/reads -type f -name '*R2.fastq.gz' | sort -V))

# Pull out a single element from the R1 and R2 arrays
CURRENT_FWD=${SAMPLES_F[${PBS_ARRAYID}]}
CURRENT_REV=${SAMPLES_R[${PBS_ARRAYID}]}

# Get the sample name for the output file
SAMPLE=$(basename ${CURRENT_FWD} | cut -f 1 -d '_')

# Run the bwa command
bwa mem \
    -t 8 -k 12 -M \
    /home/group/shared/references/ref.fa \
    ${CURRENT_FWD} ${CURRENT_REV} \
    > /panfs/roc/scratch/konox006/aln/${SAMPLE}.sam
