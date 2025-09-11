"""
Adds INFO/VAF field to SV mosaicism caller results
"""
import sys
import argparse
import truvari

def severus_vaf(entry):
    """
    Severus already calculates it
    """
    return entry.samples[0]['VAF']

def sniffles_vaf(entry):
    """
    Calculate sniffles' predicted vaf
    """
    dv = entry.samples[0]['DV']
    dr = entry.samples[0]['DR']
    tot = dv + dr
    vaf = None if tot == 0 else dv / tot
    return vaf

def sawfish_vaf(entry):
    """
    Calculate sawfish' predicted vaf
    """
    dr, dv = entry.samples[0]['AD']
    try:
        tot = dr + dv
        return dv / (dv+dr)
    except Exception:
        return 0

def svdss_vaf(entry):
    """
    Calculate sawfish' predicted vaf
    """
    return entry.info['WEIGHT'] / entry.info['COV']

def pbsv_vaf(entry):
    return sawfish_vaf(entry)

def delly_vaf(entry):
    dv = entry.samples[0]['DV'] + entry.samples[0]['RV']
    dr = entry.samples[0]['DR'] + entry.samples[0]['RR']
    tot = dv + dr
    vaf = 0 if tot == 0 else dv / tot
    return vaf

def svaba_vaf(entry):
    ad = entry.samples[0]['AD']
    tot = entry.samples[0]['DP']
    return ad / tot

def gridss_vaf(entry):
    vf = entry.samples[0]['VF']
    ref = entry.samples[0]['REF']
    refpair = entry.samples[0]['REFPAIR']
    tot = vf + ref
    big_tot = tot + refpair
    if entry.var_size() < 1000:
        vaf = vf / tot
    else:
        vaf = vf / big_tot
    return vaf

tool_vaf = {'sniffles': sniffles_vaf,
            'svdss': svdss_vaf,
            'sawfish': sawfish_vaf,
            'severus': severus_vaf,
            'pbsv': pbsv_vaf,
            'delly': delly_vaf,
            'svaba': svaba_vaf,
            'gridss': gridss_vaf,
            }

def parse_args(args):
    """
    """
    parser = argparse.ArgumentParser(prog="add_vaf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=str,
                        help="Input VCF file")
    parser.add_argument("-p", "--program", type=str, required=True, choices=tool_vaf.keys(),
                        help="Program which produced the comparison VCF")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output VCF (%(default)s)")
    args = parser.parse_args(args)
    return args

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    vcf = truvari.VariantFile(args.vcf)
    header = vcf.header.copy()
    header.add_line(('##INFO=<ID=VAF,Number=R,Type=Float,'
                     'Description="Variant Allele Fraction">'))
    # named file option
    out = truvari.VariantFile(args.output, 'w', header=header)
    for entry in vcf:
        entry.translate(header)
        alt_vaf = tool_vaf[args.program](entry)
        entry.info['VAF'] = (1-alt_vaf, alt_vaf)
        out.write(entry)
    out.close()
