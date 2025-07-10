comp=$1
bed=$2
prog=$3
base=/Users/english/code/smaht/example_bench/smaht.mims.sv.vcf.gz
ref=/Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa
name=$(basename $comp)
name=bench_${name%.vcf.gz}

tabix -f $comp
truvari bench -p 0.9 -P 0.9 --passonly --pick multi \
    --reference $ref --includebed $bed -b $base -c $comp \
    -o ${name}/

python /Users/english/code/SMaHT_MIMS/scripts/make_mosaic_summary.py -i ${name} -p $prog -o ${name}/mosaic.summary.txt
