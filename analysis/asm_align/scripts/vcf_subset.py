import sys
import pysam
import truvari

in_vcf_fn = sys.argv[1]
sex = sys.argv[2]
in_bed_fn = "/users/u233287/scratch/smaht/smaht.regions.bed"
in_par_fn = "/users/u233287/scratch/smaht/scripts/GRCh38.chrX.PAR.bed"

if sex not in ['m', 'f']:
    sys.stderr.write("Arg2 must be m|f\n")
    exit(1)

vcf = pysam.VariantFile(in_vcf_fn)

out = pysam.VariantFile("/dev/stdout", 'w', header=vcf.header)

regions = truvari.build_region_tree(vcf, includebed=in_bed_fn)
par, _ = truvari.build_anno_tree(in_par_fn)

for variant in truvari.region_filter(vcf, regions):
    if sex == 'm' and truvari.entry_within_tree(variant, par):
        variant.samples[0]['GT'] = (int(1 in variant.samples[0]["GT"]))
    elif sex == 'f' and variant.chrom == 'chrY':
        continue
    out.write(variant)
