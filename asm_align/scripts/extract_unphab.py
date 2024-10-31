import sys
import pysam
import truvari

unphab_fn, phab_reg_fn = sys.argv[1:]

vcf = pysam.VariantFile(unphab_fn)
tree, _ = truvari.build_anno_tree(phab_reg_fn)
out = pysam.VariantFile("/dev/stdout", 'w', header=vcf.header)

for entry in truvari.region_filter(vcf, tree, inside=False):
    out.write(entry)


