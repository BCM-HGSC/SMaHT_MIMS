import truvari
import numpy
import pysam

in_vcf = "/Users/english/code/SMaHT_MIMS/benchmark_v1.1/smaht.single.nophab.vcf.gz"
in_bed = "../running_qdpi/adotto.v2.ns1nointerautosome.mimstier1.bed.gz"

vcf = truvari.VariantFile(in_vcf, 'r')

cur_region = None
m_lengths = None

names = list(vcf.header.samples)
samples = []
for i in names:
    samples.append(i + '.1')
    samples.append(i + '.2')

print("chrom\tstart\tend\t%s" % ("\t".join(samples)))
for entry, region in vcf.fetch_bed(in_bed, with_region=True):
    m_region = f"{region[0]}\t{region[1].begin}\t{region[1].end - 1}"
    if cur_region is None:
        m_lengths = numpy.zeros(len(entry.samples) *  2, dtype=int)
        cur_region = m_region
    elif cur_region != m_region:
        print(cur_region, "\t".join(list(m_lengths.astype(str))), sep='\t')
        cur_region = m_region
        m_lengths = numpy.zeros(len(entry.samples) *  2, dtype=int)
    
    var_length = entry.var_size()
    var_length *= -1 if entry.var_type() == truvari.SV.DEL else 1
    for idx, sample in enumerate(entry.samples):
        m_gt = entry.samples[sample]['GT']
        m_lengths[2 * idx] += m_gt[0] * var_length
        m_lengths[2 * idx + 1] += m_gt[1] * var_length

print(cur_region, "\t".join(list(m_lengths.astype(str))), sep='\t')
