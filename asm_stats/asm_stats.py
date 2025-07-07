import sys

def read_fasta(stdin):
    name = None
    cur_len = 0
    seq_lines = []
    for line in stdin:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                yield name, cur_len
            name = line[1:].split()[0]
            cur_len = 0
        else:
            cur_len += len(line)
    if name:
        yield name, cur_len

def compute_stats(lengths):
    total = sum(lengths)
    sorted_lengths = sorted(lengths, reverse=True)

    n50 = 0
    l50 = 0
    cumulative = 0
    for i, l in enumerate(sorted_lengths):
        cumulative += l
        if cumulative >= total / 2:
            n50 = l
            l50 = i + 1
            break

    return {
        'num_sequences': len(lengths),
        'total_length': total,
        'min_length': min(lengths),
        'max_length': max(lengths),
        'n50': n50,
        'l50': l50
    }

def main():
    lengths = []

    for _, seqlen in read_fasta(sys.stdin):
        lengths.append(seqlen)

    stats = compute_stats(lengths)

    print(f"Number of sequences: {stats['num_sequences']}")
    print(f"Total length:        {stats['total_length']}")
    print(f"Minimum length:      {stats['min_length']}")
    print(f"Maximum length:      {stats['max_length']}")
    print(f"N50:                 {stats['n50']}")
    print(f"L50:                 {stats['l50']}")

if __name__ == "__main__":
    main()



