import sys
import pysam
import truvari
import itertools
import numpy as np


vcf_fn, bed_fn = sys.argv[1:]
vcf = pysam.VariantFile(vcf_fn)
out_vcf = pysam.VariantFile("/dev/stdout", 'w', header=vcf.header)

region_tree = truvari.build_region_tree(vcf, includebed=bed_fn)

vcf_i = truvari.region_filter(vcf, region_tree)

matcher = truvari.Matcher()
matcher.params.sizemin = 0
matcher.params.sizefilt = 0
matcher.params.sizemax = 100_000
chunks = truvari.chunker(matcher, ('base', vcf_i))

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

