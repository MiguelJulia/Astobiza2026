#!/bin/bash

## Unsupported java version installed by conda can cause errors, use the one in default
# conda activate base

## Download fasta and gtf
VERSION=108
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz

## Unzip references because cellranger can't index gz
gunzip Homo_sapiens.GRCh38.108.gtf.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

## Run the pipeline
nextflow run nf-core/scrnaseq --input sample_sheet.csv --genome GRCh38 --genome_fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --gtf Homo_sapiens.GRCh38.108.gtf --protocol 10XV2 --aligner cellranger -profile docker --outdir analysis --max_memory 58GB -resume

