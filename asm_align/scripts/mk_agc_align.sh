mkdir -p initial_alignments
# Hard coded :(
AGCFILE=/users/u233287/scratch/smaht/HPRC-yr1.agc
REF=/users/u233287/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p jobs
NAMES=$(~/scratch/misc_software/agc-2.1_x64-linux/agc listset $AGCFILE | grep -w "HG00438\|HG002\|HG02257\|HG02486\|HG02622\|HG005")
for i in $NAMES
do 
    echo "#!/bin/bash" > jobs/aln_hprc_${i}.sh
    bash $DIR/agc_map_haplo.sh \
        $AGCFILE \
        $i \
        $REF \
        initial_alignments/ >> jobs/aln_hprc_${i}.sh
done
