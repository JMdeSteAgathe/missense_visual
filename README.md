# Gene Variant Pathogenicity Dashboard

This Dash app queries the Ensembl REST API for a variant,
parses gene and VCF data, and visualizes missense variant scores
from gnomAD and ClinVar.

Currently, only `REVEL`, `AlphaMissense` and `CADDv1.7` are available.

As the `MPC2` is not available in ensembl [VEP api](https://rest.ensembl.org/documentation/info/vep_hgvs_post), 
it is only displayed for gnomad and clinvar variants.

## Setup

- Download the clinvar vcf.gz and vcf.gz.tbi
- Download the gnomad vcf.gz (1.3G) drive link here: [https://1drv.ms/u/c/0c297c45bbbdee7e/IQBZ_wNfIyQjTJJAQxwkFC4pAQ-mOmjK4qYSMmZ8jKdmrLk?e=5jMceC](https://1drv.ms/u/c/0c297c45bbbdee7e/IQBZ_wNfIyQjTJJAQxwkFC4pAQ-mOmjK4qYSMmZ8jKdmrLk?e=5jMceC)
- Download the gnomad vcf.gz.tbi
- Download the app.py script

## usage
