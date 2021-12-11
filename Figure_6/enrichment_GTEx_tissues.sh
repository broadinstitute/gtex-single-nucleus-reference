#!/bin/bash

## Example bash script used to perform ECLIPSER enrichment analysis for each cell type, trait and tissue combination for 8 GTEx snRNA-seq tissues and 21 complex traits (Eraslan et al, BioRxiv 2021).

## The ECLIPSER scripts can be found at: https://github.com/segrelabgenomics/ECLIPSER/
# See README.md for required python libraries, e.g.:
# pip install statsmodels

# The GTEx_clumping_files are the outputs from the clumping.py code and contain the LD clumped GWAS loci per complex trait

## Run jobs sequentially
# The prefix of the '*consolidated_clumping.tsv' files is the trait name as it appears in the GTEx_traits/ folder
for clumping_file in /example/GTEx_clumping_files/*/*consolidated_clumping.tsv 
do
for tissue_file in /example/GTEx_clumping_files/*
do  

    mkdir ../out_ECLIPSER_GTEx/

    folder=$(echo $tissue_file |  rev | cut -f1- -d '/' | rev)  
    trait=$(echo $clumping_file | sed 's/ /_/g' | rev | cut -f1 -d'/' | rev | sed 's/_consolidated_clumping.tsv//g')
    tissue=$(echo $folder | rev | cut -d '/' -f1 | rev)

    test=$(echo $clumping_file | rev | cut -f1 -d'/' | rev)
    if [[ "$test" != "Null_consolidated_clumping.tsv" ]]
    then
    python -u /src/enrichment.py \
      --sc_path /data/GTEx_snRNAseq_DGE_broad_cell_types.csv \
      --clumps_file $clumping_file  \
      --background_file $folder/Null_consolidated_clumping.tsv \
      --out_directory ../out_ECLIPSER_GTEx/ \
      --trait $trait \
      --tissue $tissue \
      --method Bayesian
    fi 
done 
done


## Run jobs in parallel using LSF job scheduler
for clumping_file in /example/GTEx_clumping_files/*/*consolidated_clumping.tsv 
do
for tissue_file in /example/GTEx_clumping_files/*
do  

    mkdir ../out_ECLIPSER_GTEx/

    folder=$(echo $tissue_file |  rev | cut -f1- -d '/' | rev)  
    trait=$(echo $clumping_file | sed 's/ /_/g' | rev | cut -f1 -d'/' | rev | sed 's/_consolidated_clumping.tsv//g')
    tissue=$(echo $folder | rev | cut -d '/' -f1 | rev)

    test=$(echo $clumping_file | rev | cut -f1 -d'/' | rev)
    if [[ "$test" != "Null_consolidated_clumping.tsv" ]]
    then
    bsub -J test -q short -n 1 -M 1000 -R rusage[mem=1000] -e err -o out "python -u /src/enrichment.py \
      --sc_path /data/GTEx_snRNAseq_DGE_broad_cell_types.csv \
      --clumps_file $clumping_file  \
      --background_file $folder/Null_consolidated_clumping.tsv \
      --out_directory ../out_ECLIPSER_GTEx/ \
      --trait $trait \
      --tissue $tissue \
      --method Bayesian"
    fi 
done 
done
