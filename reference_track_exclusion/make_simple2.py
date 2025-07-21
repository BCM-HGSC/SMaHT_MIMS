"""
Identify regions where there's ≥2 SVs over 98% similar
"""
import sys
import pysam
import truvari
import itertools
import numpy as np


vcf_fn, bed_fn = sys.argv[1:]
#out_vcf = pysam.VariantFile("/dev/stdout", 'w', header=vcf.header)
params = truvari.VariantParams(sizemin=50, 
                               sizefilt=50,
                               pctsize=0.98,
                               pctseq=0.98,
                               refdist=250,
                               sizemax=100_100)
vcf = truvari.VariantFile(vcf_fn, params=params)
region_tree = truvari.build_region_tree(vcf, includebed=bed_fn)

vcf_i = vcf.fetch_regions(region_tree)

chunks = truvari.chunker(params, ('base', vcf_i))

for chunk, _ in chunks:
    all_gts = []
    pos = []
    passes = True
    for i in range(len(chunk['base']) - 1):
        for j in range(i, len(chunk['base'])):
            a_entry = chunk['base'][i]
            b_entry = chunk['base'][j]
            if a_entry.match(b_entry).state:
                passes = False
            pos.append(a_entry.start)
            pos.append(a_entry.stop)
            pos.append(b_entry.start)
            pos.append(b_entry.stop)

            if not passes:
                break
        if not passes:
            break
    if not passes:
        start = min(pos)
        end = max(pos)
        print(chunk['base'][0].chrom, start, end, sep='\t')

