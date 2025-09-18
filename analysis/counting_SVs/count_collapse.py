import sys
import truvari
vcf_fn, bed_fn = sys.argv[1:]
vcf = truvari.VariantFile(sys.argv[1])
cnt = 0
for entry in vcf.fetch_bed(bed_fn):
    sz = entry.var_size()
    if sz < 50 or sz >= 50000:
        continue
    col = entry.info['NumCollapsed']
    if col:
        cnt += col + 1
print(cnt)
