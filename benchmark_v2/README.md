SMaHT MIMS SV Benchmark v2.0
============================

Usage
-----

Given a baseline VCF and BED file, run truvari on your `comp` result with Truvari v5.3:

```
truvari bench -f ${reference} -b ${baseline} --includebed ${bed} -c ${comp} -o result/
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

Choosing parameters
-------------------
For most callers, the default `--pctseq`, `--pctsize`, and `--refdist`  are reasonable. By providing a reference, symbolic alleles
such as `<DEL>` are automatically 'filled in' with their sequences, so there's no need to turn down sequence similarity.

For `--pick`, the default of `single` MUST be used as `multi` will artificially inflate recall (e.g. 1 comparison call matching 
to multiple baseline calls), and `ac` is only relevant for germline callers with a diploid GT.

Optionally, `--dup-to-ins` can be used to attempt to match comparison callsets' `<DUP>` into `<INS>`. By default,
Truvari will compare BNDs to symbolic and resolved SVs. BNDs can optionally be excluded with `--bnddist -1`.
