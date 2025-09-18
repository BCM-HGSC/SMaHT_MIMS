comp=$1
base=$2
bed=$3
prog=$4
ref=/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa
name=$(basename $comp)
name=new_bench_${name%.vcf.gz}

tabix -f $comp
truvari bench -p 0 -P 0.7 --passonly --pick single \
    --dup-to-ins \
    --reference $ref --includebed $bed -b $base -c $comp \
    -o ${name}/

python /Users/english/code/SMaHT_MIMS/scripts/surbscore.py --preset mims ${name} -o ${name}/score.tsv -c ${name}/counts.tsv
#truvari vcf2df -i -f -b ${name}/ ${name}/data.jl
# laytr command for making the final table?
# No, I'll deal with that in a moment... I think
