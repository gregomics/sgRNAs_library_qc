# sgRNA library diversity QC

sgRNA library were produced using Brie library from [https://www.addgene.org/pooled-library/broadgpp-mouse-knockout-brie/](addgene).

There is a total of 78,637 targeting 19,674 genes around 3+ sgRNA per gene. 

The goal of this analysis is to define if the sgRNA library cover the whole transcriptome.
4 different samples were sequenced. They should be considered as technical reps.
Here is a detailed info for the Brie library: [library content](https://media.addgene.org/cms/filer_public/be/cd/becdf7c4-ea7a-41a3-96c9-ef2bc3c85979/broadgpp-brie-library-contents.txt).
## Sequencing 

The sequencing was performed at Auckland on NextSeq 550. We received 2 zips:
  
  * FASTQ_Generation_2019-01-18_07_06_32Z-151450825.zip
  * FASTQ_Generation_2019-01-23_10_49_14Z-153142423.zip

Samples were mixed and samples on 4 lanes on 2 flowcells.


## Demultiplexing

The samples will have a mixture of i5 indexes with variable length stagger sequences (to increase template diversity). The i5 indexes will overlap between samples, thus demultiplexing uses unique i7 indexes only. This should happen on-instrument, so the only action required is to concatenate the four fastq files for each lane for each sample using the cat function. Only the R1 is required; R2 can be discarded. It's good to run fastqc at this stage of course.

Note that fastqc will show (lot of) errors since due to  low diversity library (amplicon).

## Merging per sample:

We don't need to analyse lane by lane but the whole pool.

This command will concatenated all fastq from the same sample for R1:

```
find /mnt/hcs/dsm-pathology-hibma-research-lab/crispr_screen/data/ -type f -name "*_R1_001.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,2  ; done | sort | uniq | while read P; do echo "processing Prefix: $P"; find ./ -type f -name "${P}_*_R1_001.fastq.gz"  -exec cat '{}' ';'  > ${P}.merged_R1.fq.gz ; done

```

Processing R2 as well (although we will not use them):

```
find /mnt/hcs/dsm-pathology-hibma-research-lab/crispr_screen/data/ -type f -name "*_R2_001.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,2  ; done | sort | uniq | while read P; do echo "processing Prefix: $P"; find ./ -type f -name "${P}_*_R2_001.fastq.gz"  -exec cat '{}' ';'  > ${P}.merged_R2.fq.gz ; done

```


## Trimming

sgRNA are integrated in a specific site where we know the flanking sequences.
We have done this using cutadapt and with zero mismatch tolerance.
The sequences to be trimmed are the following:

  * 5' CTTGTGGAAAGGACGAAACACCG 
  * 3' GTTTTAGAGCTAGAAATAGCAAG 

we will use cutadapt with the following option:
  
  * -g GCTTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAAATAGCAA

    This tells cutadapt that there is a linked adaptor and that we need to extract the sequence in between.

  * --trimmed-only: keep only the trimmed
  * --error-rate=0: No error rate allowed
  * -m 20 and -M 20 min and max length 20.

Here is the slurm script to launch:

```
sbatch --array=1-4%4 /mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/scripts/trimming_sgRNA_fastq.sl

```
The output of this script is stored: /mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/trimmedLinkedBC/
2 files by sample:

  * sampleName.merged_trimmed_R1.fq.gz: trimmed R1
  * sampleName.merged_trimmed_R1_only_20.fq: trimmed R1 with a remaining size of 20nts exactly.

## Counting the number of read after trimming:

Get a csv file with sample name and total number of sgRNA:

```
# creating the summary file:
touch summary_sgRNAs.csv

for fastq in /mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/trimmedLinkedBC/merged_trimmed_R1_only_20.fq
do
   # get only the filename
   fastqBase=$(basename $fastq)
   # get rid of the .merged_trimmed_R1_only_20.fq
   sample=${fastqBase/.merged_trimmed_R1_only_20.fq/}
   # counting the number of line and /4 -> 4 lines for a seq.
   count=$(echo $(cat ${fastq} | wc -l)/4|bc)
   # printing in a file
   echo "$sample,$count" >> summary_sgRNAs.csv
done


```

## constructing the sgRNA library for mapping:

We need to download the data. There are 2 sets of data available:

  1. sgRNA library targeting mouse genes:

```
wget https://media.addgene.org/cms/filer_public/be/cd/becdf7c4-ea7a-41a3-96c9-ef2bc3c85979/broadgpp-brie-library-contents.txt

```
Note: This is a CR delimitated file so we need to convert to fasta file for indexing:

```
sed -e "s/\r/\n/g" broadgpp-brie-library-contents.txt | cut -f 2,7 | awk 'BEGIN{OFS="\n"}{header=">sgRNA_" $1 "_" NR; print header, $2}' > brie_sgRNAs_library.fa

```
  2. Control sgRNAs: 

```
 wget https://www.addgene.org/static/cms/filer_public/5c/ca/5cca8516-45a7-4d83-bfb2-d960c4ab9de5/broadgpp-brie-library-controls.csv

```
convertion from csv to fasta file:


```
cat broadgpp-brie-library-controls.csv | sed -e "s/\"//g" | awk 'BEGIN{FS=","}{OFS="\n"}{header=">" $2; print header, $1}' > brie_ctrl_sgRNAs_library.fa

```

The whole sgRNA library will include both sgRNA targeting genes and control.

```
cat brie_sgRNAs_library.fa brie_ctrl_sgRNAs_library.fa > total_sgRNAs_library.fa
```

Creating the index for BWA alignment:

```
module load BWA
bwa index -a bwtsw /mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/data/total_sgRNAs_library.fa

```

## mapping on Alignment

Then align in bwa using the backtrack method and allowing nil mismatch


#Align, parse to sam, convert to bam for sorting and indexing then back to sam in case a human-readable file is needed

```
bwa_index=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/data/total_sgRNAs_library.fa
mapping_dir=/mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/mapping/
module load BWA SAMtools

for fastq in /mnt/hcs/dsm-pathology-hibma-research-lab/sgRNAs_library_qc/trimmed/*_trimmed_R1.fq.gz 
   do echo "processing fastq: $fastq"
   filename=($(basename $fastq))
   sai=$mapping_dir${filename/_trimmed_R1\.fq\.gz/\.sai} 
   bwa aln $bwa_index $file > $sai
   sam=${sai/\.sai/\.sam}
   # keep the sam for now:
   bwa samse $bwa_index $sai $fastq > $sam
   bam=${sam/\.sam/\_sorted.bam}
   samtools view -b -@ 8 $sam | samtools sort -@ 8 -o $bam
   samtools index $bam
done
```


