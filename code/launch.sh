#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# Dataset URL
DATASET_URL='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147507&format=file&file=GSE147507%5FRawReadCounts%5FHuman%2Etsv%2Egz'

# Download dataset
wget $DATASET_URL -O ./results/GEO-GSE147507.tsv.gz

# Unzip downloaded dataset
gzip -df './results/GEO-GSE147507.tsv.gz'
