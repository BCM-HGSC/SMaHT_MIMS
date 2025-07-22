Main Scripts with potential reuse across projects beyond SMaHT MIMS.

surbscore.py
------------

This tool automates the work of investigating a callers performance against different variant feature stratifications.
Note this tool requires `python3 -m pip install truvari>=5.3 scipy==1.10.1  statsmodels==0.13.5`
The `--output` report has columns:

| Column   | Definition                                                         |
| -------- | ------------------------------------------------------------------ |
| feature  | Stratification feature (e.g. `svtype`)                             |
| value    | Stratification value (e.g. `DEL`)                                  |
| pval     | Raw permutation test p-value                                       |
| acc      | Accuracy of calls in this stratification (e.g. `TP / (TP+FN)`)     |
| delta    | Observed mean ranking minus permutation tests' mean ranking        |
| q1       | 1% Quantile of permutation tests                                   |
| q50      | 50% Quantile of permutation tests                                  |
| q99      | 99% Quantile of permutation tests                                  |
| obs      | Number of variants observed in this stratification                 |
| adj_pval | Benjamini/Hochberg multiple test correction adjusted pvalue        |
| reject   | True or False based on rejection of null hypothesis at given alpha |

For the SMaHT MIMS SV benchmark, run on the truvari bench result with the command 
```python
python surbscore.py --preset mims bench_result_dir/
```

### How to interpret surbscore output
As an example, consider the below, truncated surbscore report

| feature	| value           | pval   | acc    | delta     | q1       | q50     | q99     | obs   | adj_pval | reject |
| --------- | --------------- | ------ | ------ | --------- | -------- | ------- | ------- | ----- | -------- | ------ |
| VAF_bin   | SomaticLow(<5%) | 0.0    | 0.140  | -111.7586 | -22.653  | -0.0418 | 22.6734 | 13184 | 0.0      | True   |
| isolated  | False           | 0.0    | 0.391  | -49.437   | -20.6976 | 0.0522  | 21.1852 | 12401 | 0.0      | True   |
| svtype    | INS             | 0.0182 | 0.4985 | -18.8122  | -21.0539 | -0.0179 | 21.2513 | 20545 | 0.1213   | False  |
| szbin     | [200,300)       | 0.5651 | 0.5299 | 2.3111    | -31.8288 | -0.0405 | 31.6662 | 2546  | 1.0      | False  |
| svtype    | DEL             | 0.9806 | 0.6264 | 18.8122   | -21.1796 | -0.0359 | 21.0722 | 13595 | 1.0      | False  |
| VAF_bin   | GermlineHom     | 1.0    | 0.9446 | 72.259    | -24.9471 | -0.0244 | 24.533  | 5308  | 1.0      | False  |

This report tells us that from top-to-bottom (lowest to highest p-value), we have ranked features that are most
dependent on the accuracy of the caller. For example, the top hit of `VAF_bin SomaticLow(<5%)` indicates that the
average accuracy of variants with this feature is 0.140, which is significantly lower than the average performance of
the caller across other feature/value stratifications. Furthermore, because the multiple-test correction adjusted
p-value (`adj_pval`) is still below our alpha (0.01), we can reject the null hypothesis that the caller's performance is
independent of this stratification.

While few stratifications will reach the level of significance, we can still use the sorted p-values to roughly assess a
callers performance on insignificant stratifications. For example, the third lowest p-value (`svtype INS, p=0.1213`) has
an `acc=0.4985`, which is lower than other stratifications with higher p-values (e.g. `svtype DEL, p=0.9806, acc=0.62`).
This ordering tells is that this example caller performs better on DEL than INS.

Advanced options are available in `surbscore.py --help`.

dynamic_collapse.py
-------------------

Performs SV collapsing (merging) using a dynamic threshold, updating SV genotypes. For the MIMS project, this was used
to create the 'easy' benchmark via:
```python
python dynamic_collapse.py smaht.sv.v2.hard.min45.vcf.gz \
    | python ../asm_align/scripts/vaf_annotator.py \
    | bcftools sort \
    | truvari anno numneigh \
    | bcftools sort -O z -o smaht.sv.v2.easy.vcf.gz
```
