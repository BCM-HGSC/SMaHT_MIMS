"""
Collapse variants at the same position and of the same length ±diff
while consolidating genotypes as needed
"""
import sys
import truvari
from truvari.collapse import sort_length
from functools import cmp_to_key

SORTER = cmp_to_key(sort_length)
MAXDIFF = 0
def split_types(entries):
    d = []
    i = []
    for entry in entries:
        if entry.var_type() == truvari.SV.DEL:
            d.append(entry)
        else:
            i.append(entry)
    return sorted(d, key=SORTER), sorted(i, key=SORTER)

def consolidate(a, b):
    for samp in a.samples:
        a_gt = a.samples[samp]['GT']
        b_gt = b.samples[samp]['GT']
        n_gt = []
        for allele_a, allele_b in zip(a_gt, b_gt):
            if 1 in (allele_a, allele_b):
                n_gt.append(1)
            elif 0 in (allele_a, allele_b):
                n_gt.append(0)
            else:
                n_gt.append(None)
        a.samples[samp]['GT'] = tuple(n_gt)
        a.samples[samp].phased = True

def collapse(entries):
    pos = 1 
    #n_con = 0
    annos = [0] * len(entries)
    while pos < len(entries):
        # longer, shorter
        a, b = entries[pos - 1], entries[pos]
        if a.var_size() - b.var_size() <= MAXDIFF:
            consolidate(a, b)
            annos[pos - 1] += 1
            del(entries[pos])
            del(annos[pos])
        else:
            # previous how many it increases
            pos += 1
    return entries, annos

def flush(entries, header, out):
    dels, inss = split_types(entries)
    dels = collapse(dels)
    inss = collapse(inss)
    tot_within = 0
    for entry, anno in zip(*dels):
        entry.translate(header)
        if anno:
            tot_within += 1 + anno
        entry.info["NumCollapsed"] = anno
        out.write(entry)
    for entry, anno in zip(*inss):
        entry.translate(header)
        if anno:
            tot_within += 1 + anno
        entry.info["NumCollapsed"] = anno
        out.write(entry)
    return tot_within

if __name__ == '__main__':
    in_fn = sys.argv[1]
    bed_fn = sys.argv[2]
    vcf = truvari.VariantFile(in_fn)
    header = vcf.header.copy()
    header.add_line(('##INFO=<ID=NumCollapsed,Number=1,Type=Integer,'
                     f'Description="Number of calls within {MAXDIFF} length collapsed into this call">'))
    out = truvari.VariantFile("/dev/null", 'w', header=header)
    prev_entries = [next(vcf)]
    tot_within = 0
    tot_entries = 0
    for entry in vcf.fetch_bed(bed_fn):
        if entry.var_size() < 50:
            continue
        tot_entries += 1
        if entry.chrom == prev_entries[-1].chrom and entry.pos == prev_entries[-1].pos:
            prev_entries.append(entry)
        else:
            tot_within += flush(prev_entries, header, out)
            prev_entries = [entry]
    tot_within += flush(prev_entries, header, out)
    sys.stderr.write(f"Total Entries: {tot_entries}\n")
    sys.stderr.write(f"Total Within Thresh: {tot_within}\n")

