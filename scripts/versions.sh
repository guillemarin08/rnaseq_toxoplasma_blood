#!/bin/bash
echo "FastQC:"
apptainer exec /containers/apptainer/fastqc-0.12.1.sif fastqc --version
