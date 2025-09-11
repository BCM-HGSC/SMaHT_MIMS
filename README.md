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

A helper script in `analysis/scripts/add_vaf.py` can add an INFO/VAF annotation to results from
some SV callers. This is not required, but will help with downstream analysis. The script can be easily extended to
calculate a different VAF for different callers, see the script for details.

Because SVs can be highly similar (e.g. 100bp germline homozygous insertion and 105bp somatic insertion at the same 
mosaic locus), it is important that all SVs are processed together to ensure that the best matching SVs are found. The
SMaHT MIMS SV benchmark contains the germline and somatic mutations in a single VCF. Therefore, your comparison set of
SVs should also contain germline and somatic mutations in a single VCF. Then, Truvari will perform the SV matching and
post-processing can be performed to find somatic-specific performance. 

Choosing a benchmark
--------------------

See all benchmark files in [benchmark_v2](benchmark_v2/README.md)
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

```bash
truvari bench -b benchmark_v2/smaht_mims_sv_v2_easy.vcf.gz
    --includebed benchmark_v2/smaht_mims_sv_v2_easy_regions.bed
    -o result_dir/
```
Calculating Somatic Recall
--------------------------
MIMS was designed for single-sample SV mosaicism benchmarking. This means we expect most tools to be operating on a
single BAM and producing a single VCF, within which the germline and somatic SVs are side-by-side. The MIMS VCFs
also contain the germline and somatic SVs side-by-side. We then let truvari find the optimal matches between
baseline/comparison SVs before post-filtering to calculate performance of the somatic stratifications.

To calculate somatic recall, one simply needs to subset the baseline SVs to those with a `INFO/VAF[1]` less than 0.25.
One way to do this is to take advantage of `truvari vcf2df` by running:

```bash
truvari vcf2df -i -f -b result_dir/ result_dir/data.jl
```

This will turn your benchmarking result into a pandas DataFrame. Then, a simple python script can calculate the somatic
recall.
```python
import joblib
data = joblib.load("result_dir/data.jl")
somatic_base = data[data['state'].isin(['tpbase', 'fn']) & (data['VAF_alt'] < 0.25)]
counts = somatic_base['state'].value_counts()
print("TP:" counts['tp'])
print("FN:", counts['fn'])
print("Precision:", counts['tpbase'] / len(somatic_comp))
```

To calculate somatic precision, we need to rely on the comparison SVs' VAF annotations. Again, the script 
`analysis/scripts/make_mosaic_summary.py` can add `INFO/VAF` annotations. Assuming they're inside the comparison VCF, we
can calculate somatic precision as:

```python
somatic_comp = data[data['state'].isin(['tp', 'fp']) & (data['VAF_alt'] < 0.25)]
counts = somatic_comp['state'].value_counts()
print("TP:" counts['tp'])
print("FP:", counts['fp'])
print("Precision:", counts['tp'] / len(somatic_comp))
```

Automatic Stratification Performance
------------------------------------

Annotations inside the MIMS VCFs allow one to run Truvari's StratP test. See [documentation](https://github.com/ACEnglish/truvari/wiki/Stratp-Test)
for details. Generate a stratification report of a result's recall performance:
```bash
truvari stratp --preset mims result_dir/
```
