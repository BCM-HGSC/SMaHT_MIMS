
Workspace for SMaHT MIMS benchmark creation and analysis tools

See `benchmark_v2` for the actual benchmark. 

<a href="https://raw.githack.com/BCM-HGSC/SMaHT_MIMS/main/coverage_calculator/vaf_coverage_predictor.html">Coverage Calculator</a>


Terminology
==========


| Variant Type             | Expected VAF Range | Notes                                             |
| ------------------------ | ------------------ | ------------------------------------------------- |
| Germline homozygous      |  ~100%             | Present in all cells                              |
| Germline heterozygous    |  ~50%              | Present in all cells, diploid                     |
| Somatic (clonal)         |  ~30–50%           | Tumor or clonal somatic in bulk tissue            |
| Somatic (subclonal)      |  ~5–30%            | Present in a subset of cells                      |
| Mosaic (constitutional)  |  ~1–30%            | Depends on developmental timing                   |
| Ultra-low VAF somatic    | <5%                | Needs deep/targeted sequencing                    |

| Term           | When Mutation Occurs        | How Many Cells Affected       | VAF in Bulk WGS                    |
| -------------- | --------------------------- | ----------------------------- | ---------------------------------- |
| Germline       | Before fertilization        | All                           | ~50% or 100%                       |
| Constitutional | Shortly after fertilization | Subset across tissues         | ~1–30%                             |
| Clonal         | Early in clonal expansion   | Most of the cells in a sample | ~30–50% (or higher with copy gain) |
| Subclonal      | Later during expansion      | Minority of cells             | ~1–30%                             |
