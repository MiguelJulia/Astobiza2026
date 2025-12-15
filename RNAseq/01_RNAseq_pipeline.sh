#!/bin/bash

## Requieres nextflow and singularity installed

## Run the pipeline
nextflow run nf-core/rnaseq --input sample_sheet_2.csv --genome GRCm38 -profile singularity --outdir analysis_2 --max_memory 58GB -resume
