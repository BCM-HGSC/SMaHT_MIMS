#bed=/Users/english/code/smaht/example_bench/smaht.cleaned.regions.bed
#bed=/Users/english/code/smaht/example_bench/smaht.regions.bed
bed=/Users/english/code/smaht/example_bench/smaht.cleaned.nosplit.regions.bed
#for prog in delly gridss pbsv sawfish severus sniffles svaba svdss
for prog in sniffles
do
    cd $prog
    rm -rf bench_*
    for i in *.vcf.gz; 
    do
        if [[ "$i" == *"_ill_"* ]]; then
            bash ../run_ill_truvari.sh $i $bed $prog;
        else
            bash ../run_truvari.sh $i $bed $prog;
        fi
    done 
    cd ../
done

pfx=split
for i in fn tp-base; do truvari consistency -o ${pfx}_consis.${i}.table.txt --json */*/${i}.vcf.gz > ${pfx}_consis.${i}.summary.json; done

find * -type f -name "summary.json" | tar -cvf ${pfx}_summaries.tar -T -
