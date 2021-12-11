#!/bin/bash

## Example bash script used to perform LD clumping of GWAS variants taken from Open Targets Genetics for 21 complex diseases and traits relevant 
## to 8 snRNA-seq GTEx tissues

# The GWASvar2gene input files used for the clumping of the trait GWAS variants and gene to variant mapping can be downloaded from Zenodo and saved in the data/ folder.

# Run jobs sequentially
# The /data/GTEx_traits/ folder should contain a list of files with the prefix being the tissue name that contains all traits relevant to that tissue 
for tissue in /data/GTEx_traits/*traits.txt
do
 s="${tissue##*/}"
 s="${s%_traits.txt}"
 
 out_path="/example/GTEx_clumping_files/$s"
 
 cd /example/GTEx_clumping_files
 mkdir "$s"
 
 python -u /src/clumping.py \
   --proxy_table /data/GWASVar2gene_open_targets_proxy_table.tsv \
   --traits_table /data/GWASVar2gene_open_targets_trait_table.tsv \
   --sQTL_table /data/GWASVar2gene_open_targets_sQTL_table.tsv \
   --eQTL_table /data/GWASVar2gene_open_targets_eQTL_table.tsv \
   --opentarget /data/opentargets_variant_and_gene_table.csv \
   --significant_traits $tissue \
   -ancestry EUR \
   -out_directory $out_path 
done


# Run jobs in parallel using LSF job scheduler 
for tissue in /data/GTEx_traits/*traits.txt
do
 s="${tissue##*/}"
 s="${s%_traits.txt}"
 
 out_path="/example/GTEx_clumping_files/$s"

 cd /example/GTEx_clumping_files
 mkdir "$s"
 
bsub -J test -q big -n 1 -M 60000 -R rusage[mem=60000] -e err -o out "python -u /src/clumping.py \

   --proxy_table /data/GWASVar2gene_open_targets_proxy_table.tsv \
   --traits_table /data/GWASVar2gene_open_targets_trait_table.tsv \
   --sQTL_table /data/GWASVar2gene_open_targets_sQTL_table.tsv \
   --eQTL_table /data/GWASVar2gene_open_targets_eQTL_table.tsv \
   --opentarget /data/opentargets_variant_and_gene_table.csv \
   --significant_traits $tissue \
   -ancestry EUR \
   -out_directory $out_path"
done

