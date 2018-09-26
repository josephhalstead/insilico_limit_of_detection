## Insilico Limit of Detection

## Introduction

A prototype application for downsampling a bam until we can no longer detect known variants. This allows us to know at what coverage variants can still be discovered.

The software downsamples a bam at a series of user specified rates. It will repeat this a number of user specified times to account for the randomness of the downsampling process.

The user can the see if a particular known variant can be detected at lower depths.

At present the Mutect2 variant caller is used. Users can edit the call_variants rule within the Snakefile to change this if they desire.

## Requirements and Installation

The software is designed to run on Linux in a conda environment. Follow the below instructions to install:

First follow the instructions to install Miniconda at: https://conda.io/docs/user-guide/install/linux.html

Then type the following into your terminal at a selected installation directory:

`conda env create -f envs/main.yaml`

This should install all the requirements

## Running

There are two stages to the software:

1) Run downsampling and variant calling pipeline
2) Analyse data

The first stage consists of a snakemake pipeline. The config file 'config.yaml' should be edited to give the location of your BAM files as well as entering other variables such as the replication rate and downsampling rates to try.

The first stage can be run by typing:

`snakemake final.txt`

The second stage consists of a jupyter notebook which analyses the data. Within your terminal type:

`jupyter notebook`

You can then navigate to the notebook in utils/analysis.ipynb and follow the instructions there.

## Limitations

* Does not take into account allelle frequency.



