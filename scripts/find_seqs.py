# Made by Igor Trujnara (github.com/igik20/)
"""
This script finds subsequence ranges for a feature in a Phobius short output file.
Usage: python3 find_seqs.py [input filepath] [output filepath] [feature symbol]
Input format: Phobius short prediction 
Output format:
[sequence name]    [subsequence 1 start] [subsequence 1 end]    [...other subsequences]
[...other sequences]
"""

from phobius_field_retriever import phobius_short_pred_field_selector
import sys

def find_sequences(infile, outfile, feature):
    with open(infile, 'r') as inf:
        with open(outfile, 'w') as of:
            for line in inf:
                pred = phobius_short_pred_field_selector(line, feature, False)
                of.write(str(pred[0]) + '\t' + '\t'.join([f"{t[0]} {t[1]}" for t in pred[1]]) + '\n')

if __name__ == "__main__":
    try:
        infile, outfile, feature = sys.argv[1:4]
        find_sequences(infile, outfile, feature)
    except ValueError:
        raise SystemExit("Check arguments. Usage: python3 find_seqs.py [input filepath] [output filepath] [feature symbol]")