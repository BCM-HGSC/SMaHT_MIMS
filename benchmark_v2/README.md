SMaHT MIMS SV Benchmark v2.0
============================

Usage
-----

Given a baseline VCF and BED file, run truvari on your `comp` result with:

```
truvari bench -b ${baseline} --includebed ${bed} -c ${comp} -o result/
```

Choosing a benchmark
--------------------

Quick answer is to use `smaht_mims_sv_v2_easy.vcf.gz` and `smaht_mims_sv_v2_easy_nodenseregions.bed`.

There's two SV VCFs, hard and easy. The `hard` VCF reports every SV from the HapMap alleles. The `easy` VCF
collapses SVs which are between ~5bp-30bp similar to another SVs. The `easy` VCF is more reflective of what's
reasonable for an SV caller to report, whereas the `hard` VCF expects SV callers to be able to identify e.g. SNPs, small
indel differences between SV alleles.

Each SV has two associated bed files, `fullregions` and `nodenseregions`. The `fullregions` are GRCh38 autosomes
confidently covered by each haplotype of the HapMap samples. The `nodenseregions` excludes regions where a single
haplotype produced multiple SVs. These dense regions contain complex SVs which are difficult to match to because they're
subject to alignment/representation ambiguities.

Pair the hard/easy VCF with their appropriate full/nodense regions.
