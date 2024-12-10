"""
Make counts of 
Sample Het Hom
"""
import sys
import truvari
import pysam
from collections import defaultdict

vcf_fn, bed_fn = sys.argv[1:]


def counter(vcf_i):
    lookup = defaultdict(lambda: [0, 0])

    for entry in vcf_i:
        for sample, data in entry.samples.items():
            cnt = sum(data['GT'])
            if cnt == 1:
                lookup[sample][0] += 1
            elif cnt == 2:
                lookup[sample][1] += 1
    import json
    for key in lookup:
        print(key, *lookup[key], sep='\t')

def vaf_extract(vcf_i):
    for entry in vcf_i:
        print(entry.info['VAF'][1])

if __name__ == '__main__':
    vcf = pysam.VariantFile(vcf_fn)
    regions = truvari.build_region_tree(vcf, includebed=bed_fn)

    vcf_i = truvari.region_filter(vcf, regions)
    vaf_extract(vcf_i)
