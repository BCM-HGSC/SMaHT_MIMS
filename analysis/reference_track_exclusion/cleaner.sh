
#Given the regions:
reg=/Users/english/code/smaht/example_bench/smaht.regions.bed

# Trim out anything that's within 5kbp of a gap.
bedtools slop -b 5000 -i <(zgrep -v '_' grch38.gap.bed.gz) \
    -g /Users/english/code/regioners/../references/grch38/GRCh38_1kg_mainchrs.fa.fai \
    > gap.slop.bed

# So slop the gap, and I'll probably need the genome file as well

# And then any region which doesn't fully span one of these needs to be excluded
# So -u -f 1 of this against regions
# Then subtract
bedtools intersect -v -f 1 -a grch38.bigsimplerepeat.bed.gz -b smaht.regions.bed > toobig.bed
bedtools intersect -v -f 1 -a grch38.segdups.bed.gz -b smaht.regions.bed > tooseg.bed
# And also, grab the Germline always missing set as well
bedtools slop -b 500 -i <(cut -f1-3 ../manual_inspection/to_exclude.bed) \
    -g /Users/english/code/regioners/../references/grch38/GRCh38_1kg_mainchrs.fa.fai \
    > fn.slop.bed

cat fn.slop.bed toobig.bed tooseg.bed gap.slop.bed \
    | cut -f1-3 | bedtools sort | bedtools merge > to_exclude.bed

bedtools subtract -a $reg -b to_exclude.bed | grep -v chrX | grep -v chrY > smaht.cleaned.regions.bed

python make_simple.py ~/code/smaht/example_bench/smaht.mims.sv.vcf.gz smaht.cleaned.regions.bed > split_var_candidates.bed
bedtools slop -b 500 -i split_var_candidates.bed -g ~/code/references/grch38/GRCh38_1kg_mainchrs.fa.fai > split_var_candidates.slop.bed

bedtools subtract -a smaht.cleaned.regions.bed -b split_var_candidates.slop.bed > smaht.cleaned.nosplit.regions.bed
