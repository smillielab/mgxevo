import sys
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
from ut import seq as seq_util

def main():
    parser = argparse.ArgumentParser(description='This script filters blast results based on certain percentage and length cut off.')
    parser.add_argument('--fasta', help='Input file in FASTA format.')
    parser.add_argument('--blast', help='Blast file to be filtered.')
    parser.add_argument('--filtered_blast', help='Output file for filtered blast results.')
    args = parser.parse_args()

    qlen = {k:len(v) for k,v in seq_util.read_fst(args.fasta).items()}

    with open(args.filtered_blast, 'w') as outfile:
        for line in open(args.blast):
            line = line.rstrip().split('\t')
            qid = line[0].split()[0]
            tid = line[1].split()[0]
            pct = float(line[2])
            mlen = float(line[3])
            if pct >= 90 and mlen >= 0.9*qlen[qid]:
                outfile.write('\t'.join(line) + '\n')

if __name__ == "__main__":
    main()
