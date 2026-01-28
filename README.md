# missense-visual

Compare your variant to gnomad and clinvar before applying your PP3 criteria (or discarding it):

<img width="753" height="819" alt="image" src="https://github.com/user-attachments/assets/18dc5712-1f3c-4331-8ff2-4baf88b5e5e2" />






Currently, only `REVEL`, `AlphaMissense` and `CADDv1.7` are available. `MISTIC`, `MPC2` and `popEVE` are displayed for gnomad and clinvar, but not for custom.

As the `MPC2` is not available in ensembl [VEP api](https://rest.ensembl.org/documentation/info/vep_hgvs_post), 
it is only displayed for gnomad and clinvar variants.

## Setup

- Download the /Users/jedesain/Documents/missense/clinvar_plp_ms.fully_annotated.vcf.gz and its csi
- Download the gnomad vcf.gz (1.3G) drive link here: [https://1drv.ms/u/c/0c297c45bbbdee7e/IQDlnodDmGNTRq2kMQJ0AE0PAdXzhM5P4HHX8HtbZmxx8hc?e=TBADcU](https://1drv.ms/u/c/0c297c45bbbdee7e/IQDlnodDmGNTRq2kMQJ0AE0PAdXzhM5P4HHX8HtbZmxx8hc?e=TBADcU)
- Download the index /Users/jedesain/Documents/missense/gnomad_ms.fully_annotated.vcf.gz.csi
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
