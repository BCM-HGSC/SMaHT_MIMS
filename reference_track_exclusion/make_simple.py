import sys
import pysam
import truvari
import itertools
import numpy as np


vcf_fn, bed_fn = sys.argv[1:]
#out_vcf = pysam.VariantFile("/dev/stdout", 'w', header=vcf.header)
params = truvari.VariantParams(sizemin=0, 
                               sizefilt=0,
                               sizemax=100_100)
vcf = truvari.VariantFile(vcf_fn, params=params)
region_tree = truvari.build_region_tree(vcf, includebed=bed_fn)

vcf_i = vcf.fetch_regions(region_tree)

chunks = truvari.chunker(params, ('base', vcf_i))

for chunk, _ in chunks:
    all_gts = []
    pos = []
    for entry in chunk['base']:
        all_gts.append(list(itertools.chain(*[_['GT'] for _ in entry.samples.values()])))
        pos.append(entry.start)
        pos.append(entry.stop)
    all_gts = np.array(all_gts)

    tot = np.sum(all_gts, axis=0)
    if np.max(tot) >= 2:
        start = min(pos)
        end = max(pos)
        print(chunk['base'][0].chrom, start, end, sep='\t')

