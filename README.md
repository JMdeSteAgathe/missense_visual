# missense-visual

Compare your variant to gnomad and clinvar before applying your PP3 criteria (or discarding it):

<img width="645" height="933" alt="image" src="https://github.com/user-attachments/assets/5922933c-998c-45f4-b976-744a84acd988" />





Currently, some scores (e.g. `MISTIC`, `MPC2`, `popEVE`) are displayed for gnomad and clinvar, but not for custom (not available in ensembl [VEP api](https://rest.ensembl.org/documentation/info/vep_hgvs_post)).

## Setup

Download the 2 preannotated vcf with their index
-  `clinvar_plp_ms.fully_annotated.vcf.gz` and `clinvar_plp_ms.fully_annotated.vcf.gz.csi` (available in github rep)
-  `gnomad.PASS.BCSQ.ms.fully_annotated.vcf.gz` and its index (1.3G) drive link here: [https://1drv.ms/f/c/b3a192f976546575/IgBImZb8Qas0RKhQW9K5C7ldARwgVYtpQCSjqzl_CU_zE-A?e=Vsonz4](https://1drv.ms/f/c/b3a192f976546575/IgBImZb8Qas0RKhQW9K5C7ldARwgVYtpQCSjqzl_CU_zE-A?e=Vsonz4) (the index `gnomad_ms.fully_annotated.vcf.gz.csi` is available in the github files)

Download the `gene_coord.csv.gz` (available in github files)

Download the `app.py` script

## Usage
- Adapt the 3 file paths in the `app.py` (at the top of the code)
```
VCF_FILE = "/path/to/gnomad.PASS.BCSQ.ms.fully_annotated.vcf.gz"
CLINVAR_VCF = "/path/to/clinvar_plp_ms.fully_annotated.vcf.gz"
GENE_COORD = "/path/to/gene_coord.csv.gz"
```
- Replace the gene and variant by your gene and your variant
```
TARGET_GENE = "SCN1A"
variant = "NM_001165963.4(SCN1A):c.1060G>C"
```
- Run the app.py
