import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

import ut

fna = sys.argv[1]

for record in ut.iter_seq(fna):
    contig = record[0][1:].split()[0]
    gene = contig
    beg = 1
    end = len(record[1])
    strand = '+'
    seq = record[1]
    print('\t'.join(map(str, [contig, gene, beg, end, strand, seq])))
