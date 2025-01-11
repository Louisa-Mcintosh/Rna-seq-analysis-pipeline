#!/bin/bash
set -e  # Exit on any error
exec > >(tee "pipeline.log") 2>&1  # Log everything to pipeline.log

# Directories
RAW_DIR="rawdata"               # Raw FASTQ files
TRIMMED_DIR="trimmed"           # Trimmed FASTQ files
REFERENCE_DIR="genome"          # Reference genome directory
MAPPING_DIR="mapping"           # Mapping output directory
COUNTS_DIR="feature_counts"     # Counts output directory

# Files
GTF_FILE="$REFERENCE_DIR/Homo_sapiens.GRCh38.101.gtf"
INDEX_PREFIX="$REFERENCE_DIR/genome"

# Parameters
THREADS=4  # Adjust for your system

# Create output directories
mkdir -p $TRIMMED_DIR $MAPPING_DIR $COUNTS_DIR

echo "Starting RNA-Seq Pipeline"

### Step 1: Quality Control
echo "Performing quality control on raw data..."
for file in $RAW_DIR/*.fastq.gz
do
    fastqc $file -o $RAW_DIR
done

### Step 2: Trimming
echo "Trimming raw reads for quality..."
for file in $RAW_DIR/*_1.fastq.gz
do
    base=$(basename $file _1.fastq.gz)
    trimmomatic SE -threads $THREADS \
        $RAW_DIR/${base}_1.fastq.gz \
        $TRIMMED_DIR/${base}_trimmed.fastq.gz \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

### Step 3: Alignment with Hisat2
echo "Aligning reads to the reference genome..."
for file in $TRIMMED_DIR/*.fastq.gz
do
    base=$(basename $file _trimmed.fastq.gz)
    hisat2 -p $THREADS -x $INDEX_PREFIX -U $file -S $MAPPING_DIR/${base}.sam
    samtools view -@ $THREADS -bS $MAPPING_DIR/${base}.sam | \
        samtools sort -@ $THREADS -o $MAPPING_DIR/${base}_sorted.bam
    samtools index $MAPPING_DIR/${base}_sorted.bam
done

### Step 4: Quantification with FeatureCounts
echo "Quantifying gene expression..."
for file in $MAPPING_DIR/*_sorted.bam
do
    base=$(basename $file _sorted.bam)
    featureCounts -T $THREADS -a $GTF_FILE -o $COUNTS_DIR/${base}_counts.txt $file
done

echo "RNA-Seq Pipeline Completed Successfully!"

