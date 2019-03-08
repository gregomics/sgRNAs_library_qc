#!/bin/bash

#SBATCH --job-name MapsgRNAsOnTarget
#SBATCH --mem  10000
#SBATCH --cpus-per-task  8
#SBATCH --mail-type END,FAIL
#SBATCH --output /home/REGISTRY/gimgr31p/logs/MapSgRNA_%A_%a.out
#SBATCH --error /home/REGISTRY/gimgr31p/logs/MapSgRNA_%A_%a.err
#SBATCH --time 10:00:00


trimdir=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/trimmedLinkedBC/
mapdir=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/mapping/
bwa_index=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/data/total_sgRNAs_library.fa
# declaring array to parallelise the call at SLURM level
trimmed_files=($(ls $trimdir/*_trimmed_R1_only_20.fq))

ArrayTaskID=${SLURM_ARRAY_TASK_ID}
# watch out array index starts at 0 on bash
index=ArrayTaskID-1
# getting r1  from the array:
trimmed_file=${trimmed_files[index]}
#r2file=${R1/_R1\.fq\.gz/_R2\.fq\.gz}
# no need to find the r2 since we don't use it.

echo "searching for r1: $r1file"

if [ -f $trimmed_file ]; then
   echo "Processing $trimmed_file"
   module purge
   module load BWA SAMtools
   prefix=($(basename $trimmed_file))
   sample=${prefix/_trimmed_R1_only_20.fq/}
   sample=${sample/Brie-/}
   echo "sample: $sample"
   sai=$mapdir${prefix/_trimmed_R1_only_20.fq/\.sai} 
   echo "bwa aln -t 8 $bwa_index $trimmed_file > $sai"
   srun bwa aln -t 8 $bwa_index $trimmed_file > $sai
   if [ ! -f $sai ]; then 
      echo "ERROR missing sai: $sai"
      echo "check bwa aln cmd"
      exit 2
   else
   
      sam=${sai/\.sai/\.sam}
      # keep the sam for now:
      echo "bwa samse $bwa_index -f $sam -r \"@RG\tID:$sample\tSM:$\tLB\tPL:ILLUMINA\" $sai $trimmed_file "
      srun bwa samse $bwa_index -f $sam -r "@RG\tID:$sample\tSM:$\tLB\tPL:ILLUMINA" $sai $trimmed_file 
      if [ ! -f $sam ]; then
         echo "ERROR missing sam: $sam"
         echo "check bwa samse cmd"
         exit 2
      else
         bam=${sam/\.sam/\_sorted.bam}
         echo "samtools view -b -@ 8 $sam | samtools sort -@ 8 -o $bam -"
         srun samtools view -b -@ 8 $sam | samtools sort -@ 8 -o $bam -
         echo "samtools index $bam"
         srun samtools index $bam
      fi
   fi
else
   echo "Missing r1: $r1file"
fi


