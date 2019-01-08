#!/bin/csh
#Start looping the procedures, create destination dir:

foreach folderNprefix ($fastq)
mkdir -p ${OutputDir}/${folderNprefix}

#mapping with bwa mem
BWA7 mem -t 16 ${REFGENOME}/ucsc.hg19.fasta \
${InputDir}/${folderNprefix}_R1_001.fastq.gz \
${InputDir}/${folderNprefix}_R2_001.fastq.gz | \
samtools view -Sb - > $OutputDir/${folderNprefix}/${folderNprefix}.hg19.bam && echo \
"** mapping and bam creation done **" >> $myLogBook

#sorting and indexing
GATK4 --java-options "-Xmx16G" SortSam \
-I $OutputDir/${folderNprefix}/${folderNprefix}.hg19.bam \
-O $OutputDir/${folderNprefix}/${folderNprefix}.hg19.sorted.bam  \
-SORT_ORDER coordinate \
-CREATE_INDEX TRUE && echo \
"** Sorting bam done **" >> ${myLogBook}

#markup duplicates
GATK4 --java-options "-Xmx16G" MarkDuplicates \
-I $OutputDir/${folderNprefix}/${folderNprefix}.hg19.sorted.bam \
-O $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup.bam \
-M $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_metrics.txt && echo \
"** markdup done **" >> ${myLogBook}

#add RG
GATK4 --java-options "-Xmx16G" AddOrReplaceReadGroups \
-I $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup.bam \
-O $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG.bam \
-RGID ${IDname} \
-RGLB ${RunName} \
-RGPL illumina \
-RGPU Patients \
-RGSM ${folderNprefix} \
-CREATE_INDEX TRUE && echo \
"** bam index done **" >> ${myLogBook}

#Recalibration BQSR
GATK4 --java-options "-Xmx16G" BaseRecalibrator \
-I $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG.bam \
-R ${REFGENOME}/ucsc.hg19.fasta \
-O $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_hg19.bam_recal_table \
--known-sites ${REFGENOME}/dbsnp_138.hg19.vcf && echo \
"** BQSR table done **" >> ${myLogBook}

GATK4 --java-options "-Xmx16G" ApplyBQSR \
-R ${REFGENOME}/ucsc.hg19.fasta \
-I $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG.bam  \
-O $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_recal_hg19.bam \
--bqsr-recal-file $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_hg19.bam_recal_table && echo \
"** BQSR done **" >> ${myLogBook}

#QC purpose
GATK4 --java-options "-Xmx16G" CollectInsertSizeMetrics \
-I $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_recal_hg19.bam \
-O $OutputDir/${folderNprefix}/${folderNprefix}.insert_size_metrics.txt \
-H $OutputDir/${folderNprefix}/${folderNprefix}.insert_size_histogram.pdf \
-M 0.5 && echo \
"** InsertSize Collected **" >> ${myLogBook}

# depth of coverage
# collecting Hsmetrics
GATK4 --java-options "-Xmx16G" CollectHsMetrics \
-I $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_recal_hg19.bam \
-O $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_recal_hg19_Hsmetrics.txt \
-R ${REFGENOME}/ucsc.hg19.fasta \
-BI ${IntervalListDIR}/ShiActSeqChr.bed.interval_list \
-TI ${IntervalListDIR}/ShiActSeqChr.bed.interval_list && echo \
"** HS metrics collected **" >> ${myLogBook}

#Caller
#variant calling with Mutect2, tumor file only
GATK4 --java-options "-Xmx32G" Mutect2 \
-R ${REFGENOME}/ucsc.hg19.fasta \
-L ${BEDfile} \
-I $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG_recal_hg19.bam \
-O $OutputDir/${folderNprefix}/${folderNprefix}.hg19.mutect2.vcf \
-tumor ${folderNprefix} && echo \
"** Mutect2 done **" >> ${myLogBook}

#remove unwanted bam files
rm $OutputDir/${folderNprefix}/${folderNprefix}.hg19.sorted.bam
rm $OutputDir/${folderNprefix}/${folderNprefix}.hg19.sorted.bai

rm $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup.bam
rm $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup.bai

rm $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG.bam
rm $OutputDir/${folderNprefix}/${folderNprefix}.sorted.markdup_AG.bai


#close the looping started with foreach
end
