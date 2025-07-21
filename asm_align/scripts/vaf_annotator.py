import sys
import pysam

props = {"HG00438": 0.5 / 100,
         "HG002":   2 / 100,
         "HG02257": 2 / 100,
         "HG02486": 2 / 100,
         "HG02622": 10 / 100,
         "HG005":   83.5 / 100}

def make_vaf(entry):
    """
    for each allele index (len(entry.alts) + 1), 
    """
    #occur_cnt = [0] * (len(entry.alts) + 1)
    occur_pct = [0] * (len(entry.alts) + 1)
    for samp in entry.samples:
        gt = entry.samples[samp]['GT']
        for i in gt:
            #occur_cnt[i] += 1
            occur_pct[i] += props[samp] / 2 # div2 because we're counting haploid and the proportions are diploid
    #for i in range(len(occur_cnt)):
        #occur_pct[i] = round(occur_cnt[i] / 12 * occur_pct[i], 4)
    return occur_pct

if __name__ == '__main__':
    if len(sys.argv) < 2:
        in_fn = '/dev/stdin'
    else:
        in_fn = sys.argv[1]

    vcf = pysam.VariantFile(in_fn)
    # Edit Header
    header = vcf.header.copy()
    header.add_line('##INFO=<ID=VAF,Number=R,Type=Float,'
                    'Description="Proportion of reads which should support each allele (including reference)">')
    out = pysam.VariantFile('/dev/stdout', 'w', header=header)
    for entry in vcf:
        entry.translate(header)
        entry.info['VAF'] = make_vaf(entry)
        out.write(entry)
