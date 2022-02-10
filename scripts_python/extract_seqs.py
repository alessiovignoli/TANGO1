# Made by Igor Trujnara (github.com/igik20/)
"""
This script extracts subsequences from FASTA sequences based on a tab-separated input file (format below).
Usage: python3 extract_seqs.py [range input filepath] [sequence input filepath] [output filepath]
Range input format:
[sequence name]    [subsequence 1 start] [subsequence 1 end]    [...other subsequences]
[...other sequences]
Sequence input format: canonical multiFASTA
Output format: canonical multiFASTA (1 entry per subsequence)
"""

from multiprocessing.sharedctypes import Value
import sys
import itertools

def read_seqs(seq, ints): # separated to limit indentation
    seqs = []
    current = ""
    bpts = list(itertools.chain.from_iterable(ints))
    write = False
    for i, c in enumerate(seq): # todo: what indexing in Phobius
        if (i + 1) in bpts:
            write = not write # start/stop reading sequence
            if not write: # save and reset at the end of sequence
                seqs.append(current)
                current = ""
        if write:
            current += c
    return seqs

def extract_seqs(predfile, fastafile, outfile):
    with open(fastafile, 'r') as inseqs:
        with open(predfile, 'r') as preds:
            with open(outfile, 'w') as out:
                for pred in preds:
                    l = pred.split('\t')
                    uid = l[0] # first item is the sequence ID
                    ints = [[int(a), int(b)] for [a, b] in [p.split() for p in l[1:]]] # other items are intervals
                    found = False
                    for sl in inseqs:
                        if found:
                            c = 1
                            for seq in read_seqs(sl, ints):
                                out.write(f">{uid}_{c}\n")
                                out.write(seq + '\n')
                                c += 1
                            break
                        elif sl.split()[0][1:] == uid:
                            found = True

if __name__ == "__main__":
    try:
        predfile, fastafile, outfile = sys.argv[1:4]
        extract_seqs(predfile, fastafile, outfile)
    except ValueError:
        raise SystemExit("Check arguments. Usage: python3 extract_seqs.py [range input filepath] [sequence input filepath] [output filepath]")