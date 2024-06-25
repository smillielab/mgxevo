import sys
import glob
import os
import re
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

from ut import seq as seq_util
from ut import tree as tree_util

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refmap_proc_file", help="Path to refmap processed file")
    parser.add_argument("--tree_seq_file", help="Path to tree seq file")
    parser.add_argument("--output_dir", help="Path to output directory")
    return parser.parse_args()

def create_output_dir(output_dir):
    os.makedirs(output_dir, exist_ok=True)

def read_refmap_file(refmap_proc_file):
    return pd.read_csv(refmap_proc_file, sep='\t', header=None)

def get_unique_values(df):
    return list(set(df.loc[:,1:].values.flatten().tolist()))

def get_mgseq(gi, tree_seq_file):
    gx = set(gi)
    mgseq = {}
    for record in seq_util.iter_seq(tree_seq_file):
        nm = record[0][1:].split(' ')[0]
        if nm in gx:
            mgseq[nm]= record[1].upper()
    return mgseq

def write_fasta(df, mgseq, output_dir):
    for gn in set(df.iloc[:,0]):
        dg = df.loc[df.iloc[:,0] == gn,:]
        seq = {}
        for i in dg.index:
            genes = dg.loc[i,1:]
            nm = '|'.join(['GUT_REF'] + genes.tolist())
            seq[nm] = ''.join([mgseq[gene] for gene in genes])
        tree_util.dict2fasta(seq, os.path.join(output_dir, f'{gn}.refseq.fna'))

def main():
    args = parse_arguments()
    create_output_dir(args.output_dir)
    df = read_refmap_file(args.refmap_proc_file)
    gi = get_unique_values(df)
    mgseq = get_mgseq(gi, args.tree_seq_file)
    write_fasta(df, mgseq, args.output_dir)

if __name__ == "__main__":
    main()
