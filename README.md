# sgRNA library diversity QC

sgRNA library were produced using Brie library from [https://www.addgene.org/pooled-library/broadgpp-mouse-knockout-brie/](addgene).

There is a total of 78,637 targeting 19,674 genes around 3+ sgRNA per gene. 

The goal of this analysis is to define if the sgRNA library cover the whole transcriptome.
4 different samples were sequenced. They should be considered as technical reps.
Here is a detailed info for the Brie library: [https://media.addgene.org/cms/filer_public/be/cd/becdf7c4-ea7a-41a3-96c9-ef2bc3c85979/broadgpp-brie-library-contents.txt](library content).
## Sequencing 

The sequencing was performed at Auckland on NextSeq 550. We received 2 zips:
  
  * FASTQ_Generation_2019-01-18_07_06_32Z-151450825.zip
  * FASTQ_Generation_2019-01-23_10_49_14Z-153142423.zip

Samples were mixed and samples on 4 lanes on 2 flowcells.

## Merging per sample:

This command will concatenated all fastq from the same sample for R1:

```
find /mnt/hcs/dsm-pathology-hibma-research-lab/crispr_screen/data/ -type f -name "*_R1_001.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,2  ; done | sort | uniq | while read P; do echo "processing Prefix: $P"; find ./ -type f -name "${P}_*_R1_001.fastq.gz"  -exec cat '{}' ';'  > ${P}.merged_R1.fq.gz ; done

```

Processing R2 as well:

```
find /mnt/hcs/dsm-pathology-hibma-research-lab/crispr_screen/data/ -type f -name "*_R2_001.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,2  ; done | sort | uniq | while read P; do echo "processing Prefix: $P"; find ./ -type f -name "${P}_*_R2_001.fastq.gz"  -exec cat '{}' ';'  > ${P}.merged_R2.fq.gz ; done

```

## Demultiplexing

The samples will have a mixture of i5 indexes with variable length stagger sequences (to increase template diversity). The i5 indexes will overlap between samples, thus demultiplexing uses unique i7 indexes only. This should happen on-instrument, so the only action required is to concatenate the four fastq files for each lane for each sample using the cat function. Only the R1 is required; R2 can be discarded. It's good to run fastqc at this stage of course.

Note that the fastqc will show errors since it is a low diversity library (amplicon).


## Merging sequencing lanes


## Trimming

sgRNA are integrated in a specific site where we know the flanking sequences.
We have done this using cutadapt and with zero mismatch tolerance.
The sequences to be trimmed are the following:

  * 5' CTTGTGGAAAGGACGAAACACCG (-g option)
  * 3' GTTTTAGAGCTAGAAATAGCAAG (-a option)

```
module purge
module load cutadapt/1.9.1-foss-2016b-Python-2.7.12

cutadapt -j 8 --trimmed-only --error-rate=0 -g CTTGTGGAAAGGACGAAACACCG -a GTTTTAGAGCTAGAAATAGCAAG -max-n=0 -o trimmed.fastq input.fastq

# for all the samples:

```

### Trim 5'
cutadapt --trimmed-only --error-rate=0 -g CTTGTGGAAAGGACGAAACACCG  -o output.fastq input.fastq

### Trim 3'
cutadapt --trimmed-only --error-rate=0 -a GTTTTAGAGCTAGAAATAGCAAG  -o output.fastq input.fastq

### Remove reads with N's
cutadapt --max-n=0 -o  output.fastq input.fastq

One liner to be more efficient:


#Most reads should be exactly 20 bp at this stage but an AWK script can be used to remove any that are not

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) == 20) {print header, seq, qheader, qseq}}' < input.fastq > output.fastq

## constructing the sgRNA library for mapping:

We need to download the data:

```
wget https://media.addgene.org/cms/filer_public/be/cd/becdf7c4-ea7a-41a3-96c9-ef2bc3c85979/broadgpp-brie-library-contents.txt

```

Alignment

The list of sgRNA and their corresponding sequences and target genes (also non-targeting control sgRNA) can be downloaded from Addgene:

https://www.addgene.org/pooled-library/broadgpp-mouse-knockout-brie/

This can be used to generate an index fasta for use in bwa, for example (assuming a 2-column tab delimited text file):

#In R

FileChoice<-file.choose() # choose input file
a <- read.delim(FileChoice, skip=0, sep="\t", as.is=TRUE)

fasta.out <- matrix(data=NA, nrow=2*nrow(a),ncol=1)

for(counter in 1:nrow(a)){
c.header<-(counter*2)-1
c.sequence<-(counter*2)
fasta.out[c.header,1]<-paste(">",a[counter,1],sep="")
fasta.out[c.sequence,1]<-a[counter,2]
}

write.table(fasta.out,file=paste(date(), "fastA.txt", sep=""), sep="\t", row.names=FALSE, col.names= FALSE, append=TRUE, quote = FALSE)

#In bwa

bwa index -a bwtsw /file/path/Brie_index.fa

#Then align in bwa using the backtrack method and allowing nil mismatch

#Assign a local shell variable
ind=/file/path/Brie_index.fa

#Align, parse to sam, convert to bam for sorting and indexing then back to sam in case a human-readable file is needed
bwa aln $ind file.fastq > file.sai
bwa samse $ind file.sai file.fastq > file.sam
samtools view -b -S -o file.bam file.sam
samtools sort file.bam file.sort.bam
samtools index file.sort.bam.bam
samtools view -o file.sort.sam file.sort.bam.bam

From here, it is just a matter of scoring the unique reads for each sgRNA in the index file (for example using the table function in R, but whichever implementation you favour), outputting those counts in a delimited file to work with in R or Excel (an example attached; careful not to exclude from this list any sgRNA that have no mapped reads​) and generating a histogram that shows the distribution thereof 

