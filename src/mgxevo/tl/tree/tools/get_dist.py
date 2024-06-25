import re, pickle, glob, os
import sys
import pandas as pd
import numpy as np
import itertools
from pathlib import Path
import argparse

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))

import ut
from ut import seq as seq_util
from ut import tree as tree_util

def dist(x, y):
    return np.mean([(x[i] == y[i]) or (x[i]=='N' or y[i]=='N') for i in range(len(x))])

parser = argparse.ArgumentParser()
parser.add_argument('--usher_path', help='Path to usher fasta')
parser.add_argument('--name', help='Gene id')
parser.add_argument('--output_dir', help='Output path', default='.')
args = parser.parse_args()

ifn =  f'{args.usher_path}/{args.name}.usher.fst'
fst = seq_util.read_fst(ifn)
pairs = list(itertools.combinations(list(fst.keys()),2))
dm = [dist(fst[pi[0]], fst[pi[1]]) for pi in pairs]

ii = [pi[0] for pi in pairs]
jj = [pi[1] for pi in pairs]
df = tree_util.rbind([pd.DataFrame(list(zip(ii, jj, dm))),pd.DataFrame(list(zip(jj, ii, dm)))])
df = df.pivot(index=0, columns=1, values=2)
np.fill_diagonal(df.values, 1)
df.index.name = None

ofn = f'{args.output_dir}/{args.name}.dist.tsv'
df.to_csv(ofn, sep='\t')