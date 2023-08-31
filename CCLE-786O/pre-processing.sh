#!/bin/bash
source /etc/profile
#$ -S /bin/bash
#$ -pe pvm 4

module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/samtools-1.9-gcc-5.4.0-jjq5nua
file_path=/mnt/scratch/ccle-786O
out_path=/home/ma1111/tRCC_data/trcc_wes/bam_run/
ref_path=/mnt/storage/labs/sviswanathan/GATK_ref

java -Xmx8G -jar /home/ma1111/tools/picard.jar FastqToSam \
--FASTQ $file_path/SRR8657060_1.fastq \
--FASTQ2 $file_path/SRR8657060_2.fastq \
--OUTPUT $out_path/ccle_786O_unaligned.bam \
--SAMPLE_NAME SRR8657060-786O \
--LIBRARY_NAME HC-786O_KIDNEY \
--READ_GROUP_NAME SRS4395743 \
--SORT_ORDER queryname \
--TMP_DIR $out_path/temp \

#mark adapters to mark 5' start point and add XT tag
java -Xmx8G -jar /home/ma1111/tools/picard.jar MarkIlluminaAdapters \
--INPUT $out_path/ccle_786O_unaligned.bam \
--OUTPUT $out_path/ccle_786O_markilluminaadapters.bam \
--METRICS $out_path/ccle_786O_markilluminaadapters_metrics.txt \
--TMP_DIR $out_path/temp \

#samtofastq
java -Xmx8G -jar /home/ma1111/tools/picard.jar SamToFastq \
--INPUT $out_path/ccle_786O_markilluminaadapters.bam \
--FASTQ $out_path/ccle_786O_interleaved.fastq.gz \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
--TMP_DIR $out_path/temp \

#mark adapters to mark 5' start point and add XT tag
java -Xmx8G -jar /home/ma1111/tools/picard.jar MarkIlluminaAdapters \
--INPUT $out_path/ccle_786O_unaligned.bam \
--OUTPUT $out_path/ccle_786O_markilluminaadapters.bam \
--METRICS $out_path/ccle_786O_markilluminaadapters_metrics.txt \
--TMP_DIR $out_path/temp \

#samtofastq
java -Xmx8G -jar /home/ma1111/tools/picard.jar SamToFastq \
--INPUT $out_path/ccle_786O_markilluminaadapters.bam \
--FASTQ $out_path/ccle_786O_interleaved.fastq.gz \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
--TMP_DIR $out_path/temp \

##BWA alignment
/mnt/storage/apps/BWA/bwa-0.7.17/bwa mem -K 100000000 -M -t 7 -p  $ref_path/Homo_sapiens_assembly38.fasta \
$out_path/ccle_786O_interleaved.fastq.gz | samtools view -h -b - > $out_path/ccle_786O.aligned.bam \

## Merge alignment
java -Xmx32G -jar /home/ma1111/tools/picard.jar MergeBamAlignment \
-R $ref_path/Homo_sapiens_assembly38.fasta \
--UNMAPPED_BAM $out_path/ccle_786O_unaligned.bam \
--ALIGNED_BAM $out_path/ccle_786O.aligned.bam \
--OUTPUT $out_path/ccle_786O_merged_bamalignment.bam \
--CREATE_INDEX false \
--ADD_MATE_CIGAR true \
--CLIP_ADAPTERS true \
--CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY BestMapq \
--ATTRIBUTES_TO_RETAIN XS \
--TMP_DIR $out_path/temp \

#sort sam
java -Xmx16G -jar /home/ma1111/tools/picard.jar SortSam \
-I $out_path/ccle_786O_merged_bamalignment.bam \
-O $out_path/ccle_786O_merged_bamalignment_sorted.bam \
--SORT_ORDER queryname \
--VERBOSITY INFO \

#mark duplicates
/mnt/storage/apps/gatk-4.3.0.0/gatk MarkDuplicates \
--java-options -Xmx64G \
--TMP_DIR "temp" \
--ASSUME_SORT_ORDER "queryname" \
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
--VALIDATION_STRINGENCY SILENT \
-I $out_path/ccle_786O_merged_bamalignment_sorted.bam \
-O $out_path/ccle_786O_merged_sorted_marked_duplicates.bam \
-M $out_path/ccle_786O_marked_dup_metrics.txt \
--CREATE_MD5_FILE false \

## adding RG sample name again
java -Xmx8G -jar /home/ma1111/tools/picard.jar AddOrReplaceReadGroups \
-I $out_path/ccle_786O_merged_sorted_marked_duplicates.bam \
-O $out_path/ccle_786O_rg.bam \
--RGLB HC-786O_KIDNEY \
--RGPL ILLUMINA \
--RGPU SRR8657060 \
--RGSM 786O_KIDNEY \
--CREATE_MD5_FILE true \
--CREATE_INDEX true \

## create index with samtools
/mnt/storage/apps/samtools/1.9/bin/samtools sort $out_path/ccle_786O_rg.bam -o $out_path/ccle_786O_rg_sorted.bam \
/mnt/storage/apps/samtools/1.9/bin/samtools index $out_path/ccle_786O_rg_sorted.bam 
