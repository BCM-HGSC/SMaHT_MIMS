SMaHT MIMS SV Benchmark v2.0
============================


<a href="https://raw.githack.com/BCM-HGSC/SMaHT_MIMS/main/coverage_calculator/vaf_coverage_predictor.html">Coverage Calculator</a>

Usage
-----

Given a baseline VCF and BED file, run truvari on your `comp` result with Truvari v5.3:

```
truvari bench -f ${reference} -b ${baseline} --includebed ${bed} -c ${comp} -o result/
```

Preparing your input
--------------------

The SMaHT MIMS is a benchmark for SV mosaicism. It assumes you'll be working on finding somatic mutations within e.g. a
single tissue. This means that you'll be relying on your SV caller to detect both germline and somatic mutations. This
is different from e.g. a tumor/normal pair where a subtraction pipeline can be performed to find the tumor-specific
mutations. 

Because SVs can be highly similar (e.g. 100bp germline homozygous insertion and 105bp somatic insertion at the same 
mosaic locus), it is important that all SVs are processed together to ensure that the best matching SVs are found. The
SMaHT MIMS SV benchmark contains the germline and somatic mutations in a single VCF. Therefore, your comparison set of
SVs should also contain germline and somatic mutations in a single VCF. Then, Truvari will perform the SV matching and
post-processing can be performed to find somatic-specific performance. 

Choosing a benchmark
--------------------

Quick answer is to use `smaht_mims_sv_v2_easy.vcf.gz` and `smaht_mims_sv_v2_easy_regions.bed`.

There are two SV VCFs, `hard` and `easy`. The `hard` VCF reports every SV from the HapMap alleles. The `easy` VCF
collapses SVs which are between ~5bp-30bp similar to another SVs. The `easy` VCF is more reflective of what's
reasonable for an SV caller to report, whereas the `hard` VCF expects SV callers to be able to identify e.g. SNPs, small
indel differences between SV alleles.

Each SV has two associated bed files, the default `*_regions.bed` and `complex_*_fullregions.bed`. 
The default regions are a subset of the complex fullregions where haplotypes which produce multiple nearby SVs are
excluded. Because the SMaHT MIMS SV benchmark is primarily interested in testing callers' ability to detect low-VAF
SVs, the default bed is recommended. The complex bed contains SV regions which also have low-VAF SVs, but they also
contain SVs which may not be straight forward to compare against. This means that some reported FN/FP won't truly be
false, but instead will have alignment ambiguities (e.g. split-variant representations) that prevent a 1-to-1 match from
being found, causing precision/recall to be underestimated.

Choosing parameters
-------------------
For most callers, the default `--pctseq`, `--pctsize`, and `--refdist`  are reasonable. By providing a reference, symbolic alleles
such as `<DEL>` are automatically 'filled in' with their sequences, so there's no need to turn down sequence similarity.

For `--pick`, the default of `single` MUST be used as `multi` will artificially inflate recall (e.g. 1 comparison call matching 
to multiple baseline calls), and `ac` is only relevant for germline callers with a diploid GT.

Optionally, `--dup-to-ins` can be used to attempt to match comparison callsets' `<DUP>` into `<INS>`. By default,
Truvari will compare BNDs to symbolic and resolved SVs. BNDs can optionally be excluded with `--bnddist -1`.
