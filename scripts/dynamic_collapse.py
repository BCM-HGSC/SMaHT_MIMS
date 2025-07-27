"""
Given the phab harmonized SVs from SMaHT MIMS,
collapse highly similar SVs using a dynamic threshold.
This will find SVs of the same type and start position within 100bp and highly similar size to collapse.
Variants with higher VAF are preferred when choosing a variant representation to keep.
Genotypes from HapMap haplotypes are preserved.

The size similarity threshold is determined dynamically based on the size of the smaller SV in a pair.
For example, This will allow e.g. 50bp/55bp SVs to collapse (90% similarity), and e.g. 1000bp and 
1020bp SVs to collapse (98% similarity). Without dynamic thresholding, the 50/55bp would be missed at
98% similarity. Similarly, a 1000bp would overmerge with everything ±100bp at 90% similarity. 
The exact thresholds are subject to change, but the idea remains.
"""
import sys
import truvari

def dyn_thresh(size, min_diff=5, max_diff=30, s_min=50, s_max=1500):
    if s_min is None:
        s_min = 2 * min_diff
    if s_max is None:
        s_max = 2 * max_diff

    if size <= s_min:
        allowed_diff = min_diff
    elif size >= s_max:
        allowed_diff = max_diff
    else:
        # Linear interpolation between min_diff and max_diff
        scale = (size - s_min) / (s_max - s_min)
        allowed_diff = min_diff + scale * (max_diff - min_diff)

    return 1 - allowed_diff / size

from functools import cmp_to_key
@cmp_to_key
def sort_by_vaf(b1, b2):
    """
    Order entries from highest to lowest VAF, ties are by alphanumeric of REF
    """
    s1 = b1.info['VAF'][1]
    s2 = b2.info['VAF'][1]
    if s1 < s2:
        return 1
    if s1 > s2:
        return -1
    if b1.ref < b2.ref:
        return 1
    if b1.ref > b2.ref:
        return -1
    return 0

def match(base, comp):
    """
    Collapse two calls. Update the base in-place as we go. Return True if these were collapsed
    Conditions:
    #   Types match
    #   Position within 100bp
    #   Calculate dyn_thresh based on min(size) and use that as sizesim
    #   Collapse the genotypes into the kept variant - Note we can just recalculate VAF on the other side
    #   Update the NumCollapsed
    """
    if not base.same_type(comp):
        return False
    if not abs(base.pos - comp.pos) <= 100:
        return False
    thresh = dyn_thresh(min(base.var_size(), comp.var_size()))
    if base.sizesim(comp)[0] < thresh:
        return False
    if 'NumCollapsed' in base.info:
        base.info['NumCollapsed'] += 1
    else:
        base.info['NumCollapsed'] = 1
    consolidate(base, comp)
    return True

def consolidate(a, b):
    """
    Consoliate genotypes from b into a
    """
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

if __name__ == '__main__':
    # Expecting the phab 'hard' VCF.
    # I think we might want to select >=45bp because 50bp events might be near 45 which
    # would presumably be captured by the variant. Then we can post-filter by 50bp again to remove any
    # straggler <50bp events
    params = truvari.VariantParams(sizemin=0, sizefilt=0)
    vcf = truvari.VariantFile(sys.argv[1], params=params)
    header = vcf.header.copy()
    header.add_line(('##INFO=<ID=NumCollapsed,Number=1,Type=Integer,'
                     'Description="Number of calls collapsed into this call by truvari">'))

    # Update header with NumCollapsed
    out = truvari.VariantFile("/dev/stdout", 'w', header=header)

    # Chunk it up
    chunks = truvari.chunker(params, ('base', vcf))
    for chunk, _ in chunks:
        # Prefer the higher VAF variants
        chunk['base'].sort(key=sort_by_vaf)
        while chunk['base']:
            next_call = chunk['base'].pop(0)
            next_call.translate(header)
            if next_call.var_size() < 45:
                out.write(next_call)
                continue
            idx = 0 # next call to compare to
            while idx < len(chunk['base']):
                if match(next_call, chunk['base'][idx]):
                    chunk['base'].pop(idx) # remove collapsed
                else:
                    idx += 1
            out.write(next_call)
    out.close()
    # After done, recalculate NumNeighbors
    # And then do the 'make_simple' to recreate the regions (12 will still be fine, I think)
    # With, obviously, the phabdidnthelp exclusion.
