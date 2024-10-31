# Select for confidently covered regions
# paths are relative to this script
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

TMPDIR="temp"

genome=$DIR/GRCh38_1kg_mainchrs.genome_bedtools.bed

cat sample_beds/*bed | bedtools sort | bedtools genomecov -i - -g $genome -bga > $TMPDIR/all_coverage.bed

# 6x for all chromosomes except chrY which needs 3x
awk '($4 == 6 && $1 != "chrY") || ($4 == 3 && $1 == "chrY")' $TMPDIR/all_coverage.bed > smaht.regions.bed
