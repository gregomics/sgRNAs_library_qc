#!/bin/bash

#SBATCH --job-name TrimSgRNA
#SBATCH --mem-per-cpu   10000
#SBATCH --mail-type END,FAIL
#SBATCH --output /home/REGISTRY/gimgr31p/logs/TrimSgRNA_%A_%a.out
#SBATCH --error /home/REGISTRY/gimgr31p/logs/TrimSgRNA_%A_%a.err
#SBATCH --time 10:00:00


seqdir=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/data/
outdir=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/trimmedLinkedBC/

# declaring array to parallelise the call at SLURM level
r1merged_files=($(ls $seqdir/*.merged_R1.fq.gz))

ArrayTaskID=${SLURM_ARRAY_TASK_ID}
# watch out array index starts at 0 on bash
index=ArrayTaskID-1
# getting r1  from the array:
r1file=${r1merged_files[index]}
#r2file=${R1/_R1\.fq\.gz/_R2\.fq\.gz}
# no need to find the r2 since we don't use it.

echo "searching for r1: $r1file"

if [ -f $r1file ]; then
   echo "found"
   module purge
   module load cutadapt/1.9.1-foss-2016b-Python-2.7.12
   module load Python/3.5.2-foss-2016b
   prefix=($(basename $r1file))
   prefixOut=$outdir${prefix/_R1\.fq\.gz/_trimmed_R1\.fq\.gz}
   only20nts=$outdir${prefix/_R1\.fq\.gz/_trimmed_R1_only_20.fq}
   # Note new 5': CAACTTGTGGAAAGGACGAAACACCG
   echo "running: cutadapt -j 8 --trimmed-only --error-rate=0 -g GCTTGTGGAAAGGACGAAACACCG -a GTTTTAGAGCTAGAAATAGCAAG --max-n=0 -o $prefixOut $r1file"
   
   srun cutadapt -j 8 --trimmed-only --error-rate=0 -g GCTTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAA --max-n=0 -o $prefixOut $r1file
   
   srun zcat $prefixOut | awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) == 20) {print header, seq, qheader, qseq}}' > $only20nts

   
else
   echo "Missing r1: $r1file"
fi


