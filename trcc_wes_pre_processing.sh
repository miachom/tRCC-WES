module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/samtools-1.9-gcc-5.4.0-jjq5nua

file_path=/mnt/storage/labs/sviswanathan/tRCC_data/trcc_wes/cell_line_raw_data
out_path=/home/ma1111/tRCC_data/trcc_wes/bam_run/A3
ref_path=/mnt/storage/labs/sviswanathan/GATK_ref

java -Xmx8G -jar /home/ma1111/tools/picard.jar FastqToSam \
--FASTQ $file_path/A3_CKDN230014514-1A_H735LDSX7_L1_1.fq \
--FASTQ2 $file_path/A3_CKDN230014514-1A_H735LDSX7_L1_2.fq \
--OUTPUT $out_path/A3_unaligned.bam \
--SAMPLE_NAME A3_CKDN230014514 \
--LIBRARY_NAME Lib-A3 \
--READ_GROUP_NAME sample_A3 \
--SORT_ORDER queryname \
--TMP_DIR $out_path/temp \

#mark adapters to mark 5' start point and add XT tag
java -Xmx8G -jar /home/ma1111/tools/picard.jar MarkIlluminaAdapters \
--INPUT $out_path/A3_unaligned.bam \
--OUTPUT $out_path/A3_markilluminaadapters.bam \
--METRICS $out_path/A3_markilluminaadapters_metrics.txt \
--TMP_DIR $out_path/temp \

#3samtofastq interleave
java -Xmx8G -jar /home/ma1111/tools/picard.jar SamToFastq \
--INPUT $out_path/A3_markilluminaadapters.bam \
--FASTQ $out_path/A3_interleaved.fastq.gz \
--CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true \
--TMP_DIR $out_path/temp \

## Merge alignment
java -Xmx32G -jar /home/ma1111/tools/picard.jar MergeBamAlignment \
-R $ref_path/Homo_sapiens_assembly38.fasta \
--UNMAPPED_BAM $out_path/A3_unaligned.bam \
--ALIGNED_BAM $out_path/A3.aligned.bam \
--OUTPUT $out_path/A3_merged_bamalignment.bam \
--CREATE_INDEX false \
--ADD_MATE_CIGAR true \
--CLIP_ADAPTERS true \
--CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY BestMapq \
--ATTRIBUTES_TO_RETAIN XS \
--TMP_DIR $out_path/temp \

##BWA alignment
/mnt/storage/apps/BWA/bwa-0.7.17/bwa mem -K 100000000 -M -t 7 -p  $ref_path/Homo_sapiens_assembly38.fasta \
$out_path/A3_interleaved.fastq.gz | samtools view -h -b - > $out_path/A3.aligned.bam \

#sort sam
java -Xmx16G -jar /home/ma1111/tools/picard.jar SortSam \
-I $out_path/A3_merged_bamalignment.bam \
-O $out_path/A3_merged_bamalignment_sorted.bam \
--SORT_ORDER queryname \
--VERBOSITY INFO \

#mark duplicates
/mnt/storage/apps/gatk-4.3.0.0/gatk MarkDuplicates \
--java-options -Xmx64G \
--TMP_DIR "temp" \
--ASSUME_SORT_ORDER "queryname" \
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
--VALIDATION_STRINGENCY SILENT \
-I $out_path/A3_merged_bamalignment_sorted.bam \
-O $out_path/A3_merged_sorted_marked_duplicates.bam \
-M $out_path/A3_marked_dup_metrics.txt \
--CREATE_MD5_FILE false \

## adding RG sample name again
java -Xmx8G -jar /home/ma1111/tools/picard.jar AddOrReplaceReadGroups \
-I $out_path/A3_merged_sorted_marked_duplicates.bam \
-O $out_path/A3_CKDN230014513_rg.bam \
--RGLB Lib-A1 \
--RGPL ILLUMINA \
--RGPU LDSX7 \
--RGSM A3_CKDN230014514 \
--CREATE_MD5_FILE true \
--CREATE_INDEX true \


## create index with samtools
samtools sort $out_path/A3_CKDN230014514_rg.bam -o $out_path/A3_CKDN230014514_sorted.bam \
samtools index $out_path/A3_CKDN230014514_sorted.bam \
