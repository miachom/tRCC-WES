#!/bin/bash
source /etc/profile
#$ -S /bin/bash
#$ -pe pvm 8

#module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/samtools-1.9-gcc-5.4.0-jjq5nua
file_path=/mnt/storage/labs/sviswanathan/tRCC_data/trcc_wes/cell_line_raw_data
out_path=/home/ma1111/tRCC_data/trcc_wes/bam_run/A2
ref_path=/mnt/storage/labs/sviswanathan/GATK_ref

java -Xmx8G -jar /home/ma1111/tools/picard.jar FastqToSam \
--FASTQ $file_path/A2_CKDN230014513-1A_H75FVDSX7_L3_1.fq \
--FASTQ2 $file_path/A2_CKDN230014513-1A_H75FVDSX7_L3_2.fq \
--OUTPUT $out_path/A2_unaligned.bam \
--SAMPLE_NAME A2_CKDN230014513 \
--LIBRARY_NAME Lib-A1 \
--READ_GROUP_NAME sample_A2 \
--SORT_ORDER queryname \
--TMP_DIR $out_path/temp \

#mark adapters to mark 5' start point and add XT tag
java -Xmx8G -jar /home/ma1111/tools/picard.jar MarkIlluminaAdapters \
--INPUT $out_path/A2_unaligned.bam \
--OUTPUT $out_path/A2_markilluminaadapters.bam \
--METRICS $out_path/A2_markilluminaadapters_metrics.txt \
--TMP_DIR $out_path/temp \

#samtofastq
java -Xmx8G -jar /home/ma1111/tools/picard.jar SamToFastq \
--INPUT $out_path/A2_markilluminaadapters.bam \
--FASTQ $out_path/A2_interleaved.fastq.gz \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
--TMP_DIR $out_path/temp \

##BWA alignment
/mnt/storage/apps/BWA/bwa-0.7.17/bwa mem -K 100000000 -M -t 7 -p  $ref_path/Homo_sapiens_assembly38.fasta \
$out_path/A2_interleaved.fastq.gz | samtools view -h -b - > $out_path/A2.aligned.bam \

## Merge alignment
java -Xmx16G -jar /home/ma1111/tools/picard.jar MergeBamAlignment \
-R $ref_path/Homo_sapiens_assembly38.fasta \
--UNMAPPED_BAM $out_path/A2_unaligned.bam \
--ALIGNED_BAM $out_path/A2.aligned.bam \
--OUTPUT $out_path/A2_merged_bamalignment.bam \
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
java -Xmx8G -jar /home/ma1111/tools/picard.jar SortSam \
--input $out_path/A2_merged_bamalignment.bam \
--output $out_path/A2_merged_bamalignment_sorted.bam \
--sort-order queryname \
--verbosity WARNING \

#mark duplicates
/mnt/storage/apps/gatk-4.3.0.0/gatk MarkDuplicates \
--java-options -Xmx64G \
--TMP_DIR "temp" \
--ASSUME_SORT_ORDER "queryname" \
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
--VALIDATION_STRINGENCY SILENT \
-I $out_path/A2_merged_bamalignment_sorted.bam \
-O $out_path/A2_merged_sorted_marked_duplicates.bam \
-M $out_path/A2_marked_dup_metrics.txt \
--CREATE_MD5_FILE false \

## adding RG sample name again
java -Xmx8G -jar /home/ma1111/tools/picard.jar AddOrReplaceReadGroups \
-I $out_path/A2_merged_sorted_marked_duplicates.bam \
-O $out_path/A2_CKDN230014513.bam \
--RGLB Lib-A1 \
--RGPL ILLUMINA \
--RGPU FVDSX7 \
--RGSM A2_CKDN230014513 \
--create-output-bam-md5 true \
--create-output-bam-index true \

## create index with samtools
/mnt/storage/apps/samtools/1.9/bin/samtools index $out_path/A2_CKDN230014513.bam

