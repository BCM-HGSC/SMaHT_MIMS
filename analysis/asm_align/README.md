Pipeline for mapping HPRC assemblies and creating the

1) agc_map_haplo.sh
2) covered_regions.sh and covered_regions_female.sh
3) region_intersection.sh
4) scripts/mk_mergehaps.sh to combine them
5) scripts/vcf_subset.py on the hapo_merged vcfs
6) bcftools merge the subset variants - fill missing because we know they're all covered by the alignments
bcftools merge -m none -0

7) density
python scripts/density_counter.py smaht.variants.vcf.gz > smaht.density.bed
more than 1 indel in a region

7) phab
truvari phab -r smaht.density.bed -b smaht.variants.vcf.gz  -f
~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa -t 6 --maxsize 100000 -m '--auto --thread 4'
-o smaht.variants.phab.vcf.gz

8) combine

±100bp to the smaht.density.bed because that's what phab added and it pulled in variants.
So, Then I can use `scripts/extract_unphab.py` to make the unphab vcf and then concat/sort with the phab vcfs to make
the finished product.


```
python scripts/extract_unphab.py smaht.variants.vcf.gz smaht.density.buffer.bed | bgzip > smaht.variants.unphab.vcf.gz
tabix smaht.variants.unphab.vcf.gz
bcftools concat -a smaht.variants.unphab.vcf.gz chrvcfs/*.vcf.gz | bcftools sort -O z -o smaht.singleallele.vcf.gz
```

9) Turn it into a multi-allelic VCF, and then fill-tags

```
bcftools norm -N -f ref -c s -m+any
- or -
bcftools norm -N -f ref -c s -a --atom-overlaps . -m+any
```

Trying
```
bcftools norm -T autosomes -N -f ~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa -c s -m+any
smaht.singleallele.vcf.gz  | bcftools +fill-tags -O z -o smaht.multiallelic.vcf.gz
```
Had to stick to autosomes because the haploid/diplod if chrX broke bcftools.

Also, I'm gong to go ahead and make the nonphab multiallelic vcf with
```
bcftools norm -T autosomes.bed  -N -f ~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa -c s
-m+any smaht.variants.vcf.gz  | bcftools +fill-tags -O z -o smaht.multiallelic.nophab.vcf.gz
```

10) Cleaning

Move/rename - also fix samples in header with 
````
for i in version_1.0/*.vcf.gz;
do 
	name=$(basename $i)
	bcftools reheader -h newheader.txt -o version_1.0.1/$name$i
done
```

10) Annotate VAF

I'm also taking this opportunity to clean/standardize all the other info fields
I got code somewhere that has me started. 

```
python scripts/vaf_annotator.py version_1.0.1/smaht.multiallelic.nophab.vcf.gz \
	| bcftools annotate -x FORMAT/FT,FORMAT/AD,FORMAT/BPDP,INFO/SVLEN,INFO/SVTYPE \
	| bgzip > version_1.0.2/smaht.multiallelic.nophab.vcf.gz
```


