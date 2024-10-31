import glob
from collections import defaultdict, Counter

cnts = defaultdict(Counter)
lens = defaultdict(Counter)
for i in glob.glob("initial_alignments/*/cov.bed"):
    name = i.split('/')[1][len("hprc_"):]
    with open(i) as fh:
        for line in fh:
            data = line.strip().split('\t')
            if data[3] == "0":
                state = "Uncovered"
            elif data[3] == "1":
                state = "Confident"
            else:
                state = "Multi"
            m_len = int(data[2]) - int(data[1])
            cnts[name][state] += 1
            lens[name][state] += m_len

import joblib
joblib.dump([cnts, lens], "covstats.jl")

