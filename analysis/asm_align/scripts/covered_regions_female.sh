# Select for confidently covered regions
# paths are relative to this script
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mat_cov=$1
pat_cov=$2

genome=$DIR/GRCh38_1kg_mainchrs.genome_bedtools.bed
par_reg=$DIR/GRCh38.chrX.PAR.bed

TMPDIR="temp"
OUTDIR="coverage"
mkdir -p $TMPDIR
mkdir -p $OUTDIR

# Singly covered
awk '$4 == 1' $mat_cov | cut -f1-3 > $TMPDIR/hq_mat_cov.bed
awk '$4 == 1' $pat_cov | cut -f1-3 > $TMPDIR/hq_pat_cov.bed

cat $TMPDIR/hq_mat_cov.bed $TMPDIR//hq_pat_cov.bed \
    | bedtools sort \
    | bedtools genomecov -i - -g $genome -bga > $TMPDIR/hq_diploid_cov.bed


# Everywhere requires 2x except chrY which is missing
awk '$4 == 2 && $1 != "chrY"' $TMPDIR/hq_diploid_cov.bed | cut -f1-3 > $TMPDIR/adotto_hq_cov.bed
