Run phab to harmonize
Subset to variants ≥45bp
Run this to collapse dynamically and update the annotations
python dynamic_collapse.py smaht.sv.v2.hard.min45.vcf.gz | python ../asm_align/scripts/vaf_annotator.py | bcftools sort | truvari anno numneigh | bcftools sort -O z -o collapsed.vcf.gz
