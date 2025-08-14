import sys
import truvari
vcf_fn, bed_fn = sys.argv[1:]
vcf = truvari.VariantFile(sys.argv[1])

from collections import Counter

ty_cnt = Counter()
pos_cnt = 0
prev_pos = (None, None)

for entry in vcf.fetch_bed(bed_fn):
    sz = entry.var_size()
    if sz < 50 or sz >= 50000:
        continue
    cur_pos = entry.chrom, entry.pos
    if cur_pos != prev_pos:
        pos_cnt += 1
        prev_pos = cur_pos
    ty_cnt[entry.var_type().name] += 1
import json
print("Positions:", pos_cnt)
print("Variants:", json.dumps(ty_cnt, indent=4))
