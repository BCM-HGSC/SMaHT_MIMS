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


# Non-chrX/chrY requires 1x per haplotype
awk '$4 == 2 && $1 != "chrX" && $1 != "chrY"' $TMPDIR/hq_diploid_cov.bed | cut -f1-3 > $TMPDIR/adotto_hq_cov.bed

## IF MALE
# chrY requires 1x
awk '$4 == 1 && $1 == "chrY"' $TMPDIR/hq_diploid_cov.bed | cut -f1-3 >> $TMPDIR/adotto_hq_cov.bed

# chrX requires 2x on PAR regions
# 1) take out just chrX
grep chrX $TMPDIR/hq_diploid_cov.bed > $TMPDIR/hq_diploid_cov.chrX.bed
# 2) subtract PAR from chrX coverage and get only 1x covered
bedtools subtract -a $TMPDIR/hq_diploid_cov.chrX.bed -b $par_reg \
    | awk '$4 == 1' \
    | cut -f1-3 >> $TMPDIR/adotto_hq_cov.bed

# 3) compliment PAR to make non-PAR
bedtools complement -i $par_reg -g $genome > $TMPDIR/nonpar.bed
# 4) subtract non-PAR from chrX coverage and get only 2x covered
bedtools subtract -a $TMPDIR/hq_diploid_cov.chrX.bed -b $TMPDIR/nonpar.bed \
    | awk '$4 == 2' \
    | cut -f1-3 >> $TMPDIR/adotto_hq_cov.bed
# note, these won't be sorted, but we sort in a couple steps



