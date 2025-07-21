vcf=/Users/english/code/SMaHT_MIMS/benchmark_v2/smaht.sv.v2.medium.nosmall.vcf.gz
bed=/Users/english/code/SMaHT_MIMS/benchmark_v2/smaht.v2.regions.nodense.bed
for prog in delly gridss pbsv sawfish severus sniffles svaba svdss
do
    cd $prog
    rm -rf new_bench_*
    for i in *.vcf.gz; 
    do
        if [[ "$i" == *"_ill_"* ]]; then
            bash ../run_ill_truvari.sh $i $vcf $bed $prog;
        else
            bash ../run_truvari.sh $i $vcf $bed $prog;
        fi
    done 
    cd ../
done

pfx=split
for i in fn tp-base; do truvari consistency -o ${pfx}_consis.${i}.table.txt --json */*/${i}.vcf.gz > ${pfx}_consis.${i}.summary.json; done

find * -type f -name "summary.json" | tar -cvf ${pfx}_summaries.tar -T -
