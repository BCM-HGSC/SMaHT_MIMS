comp=$1
base=$2
bed=$3
prog=$4
ref=/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa
name=$(basename $comp)
name=new_bench_${name%.vcf.gz}

tabix -f $comp
truvari bench -p 0 -P 0.7 --passonly --pick single \
    --reference $ref --includebed $bed -b $base -c $comp \
    -o ${name}/

python /Users/english/code/SMaHT_MIMS/scripts/make_mosaic_summary.py -i ${name} -p $prog -o ${name}/mosaic.summary.txt
#truvari vcf2df -i -f -b ${name}/ ${name}/data.jl
# laytr command for making the final table?
# No, I'll deal with that in a moment... I think
