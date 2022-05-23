# Made by Igor Trujnara (github.com/igik20/)
"""
This script counts the occurrences of every k-mer of set length in all subsequences in a file, grouped by sequence ID.
Usage: python3 kmer_count.py [input filepath] [output filepath] [k-mer size]
Input format: canonical multiFASTA (1 entry per subsequence)
Output format:
[sequence name]
[k-mer 1 content] [k-mer 1 count]
[...other k-mers]
"""

import sys
from collections import defaultdict as dd

def gspl(seq, n):
    return [''.join(item) for item in zip(*[seq[n:] for n in range(n)])]

def count_kmers(infile, outfile, k_len):
    res = []
    with open(infile, 'r') as inf:
        for line in inf:
            if line[0] == '>':
                l = line.rstrip().split()
                uid = l[0][1:]
            else:
                seq = line.rstrip()
                counts = dd(int)
                for kmer in gspl(seq, int(k_len)):
                    counts[kmer] += 1
                res.append( (uid, counts) )
    with open(outfile, 'w') as f:
        for it in res:
            f.write("$ " + it[0] + '\n')
            f.write(f"$ Total: {sum(it[1].values())}\n")
            f.write(f"$ Unique: {len(it[1])}\n")
            for k, v in sorted(it[1].items(), key = lambda i: i[1], reverse = True):
                f.write(f"{k} {v}\n")

if __name__ == "__main__":
    try:
        infile, outfile, k_len = sys.argv[1:4]
        count_kmers(infile, outfile, k_len)
    except ValueError:
        raise SystemExit("Check arguments. Usage: python3 kmer_count.py [input filepath] [output filepath] [k-mer size]")