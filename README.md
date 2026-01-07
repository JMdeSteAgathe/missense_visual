# missense-visual

Compare your variant to gnomad and clinvar before applying your PP3 criteria (or discarding it):

<img width="690" height="778" alt="image" src="https://github.com/user-attachments/assets/d1175cde-ff48-4469-b47f-abdd1d574994" />




Currently, only `REVEL`, `AlphaMissense` and `CADDv1.7` are available.

As the `MPC2` is not available in ensembl [VEP api](https://rest.ensembl.org/documentation/info/vep_hgvs_post), 
it is only displayed for gnomad and clinvar variants.

## Setup

- Download the clinvar vcf.gz and vcf.gz.tbi
- Download the gnomad vcf.gz (1.3G) drive link here: [https://1drv.ms/u/c/0c297c45bbbdee7e/IQBZ_wNfIyQjTJJAQxwkFC4pAQ-mOmjK4qYSMmZ8jKdmrLk?e=5jMceC](https://1drv.ms/u/c/0c297c45bbbdee7e/IQBZ_wNfIyQjTJJAQxwkFC4pAQ-mOmjK4qYSMmZ8jKdmrLk?e=5jMceC)
- Download the gnomad vcf.gz.tbi
- Download the app.py script

## Usage
Replace the gene and variant by your gene and your variant
```
TARGET_GENE = "SCN1A"
variant = "NM_001165963.4(SCN1A):c.1060G>C"
```

## Known Bugs
- hovering a clinvar or custom variant displays always empty (0) gnomad values, which is often incorrect
- 
