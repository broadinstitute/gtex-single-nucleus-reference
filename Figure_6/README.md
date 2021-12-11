## README for scripts used to run ECLIPSER and generate plots

The bash scripts in this subdirectory were used to run the two main steps of GWAS cell type specific enrichment analysis pplied to 21 complex traits by 8 GTEx snRNA-seq tissues using ECLIPSER. The ECLIPSER scripts can be found at: https://github.com/segrelabgenomics/ECLIPSER/


### Perform LD clumping of GWAS variants: `clumping_GTEx_tissues.sh`

GWAS variants that share LD proxy variants, or eGenes or sGenes were collapsed into a single locus. This was done separately for the GWAS variants for each selected traits of interest and for all null trait GWAS loci per tissue.


### Perform GWAS enrichment analysis with `enrichment_GTEx_tissues.sh`

The significance of a GWAS locus set cell type specificity score per cell type and trait is assessed against the distribution of cell type specificity scores of a background set of GWAS loci using a Bayesian Fisher’s exact test or permutation approach.
