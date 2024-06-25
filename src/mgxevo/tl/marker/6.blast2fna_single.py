import argparse, sys, glob, os, re
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import itertools

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))
from ut import seq as seq_util
from ut import tree as tree_util

def main(args):
    
    ifn = args.filtered_blast
    blast = pd.read_csv(ifn, sep='\t', header=None,low_memory=False)
    blast.columns = ['query', 'target', 'id', 'alnlen', 'mism', 'opens', 'qlo', 'qhi', 'tlo', 'thi', 'evalue', 'bits']

    refmap  = blast.dropna().groupby('query')['target'].apply(list).to_dict()

    # save refmap
    gmap = blast.groupby('query', as_index=False)['target'].apply(','.join).set_index('query').to_dict()['target']
    ids = pd.read_csv(args.multigene_table,sep='\t', usecols=['genome', *args.genes])
    ids = ids.set_index('genome')
    ids = ids.replace(gmap)
    ids.to_csv(f'{args.output_dir}/tree-seeds.mgyg.filtered.ref_map.txt', sep='\t')
    
    # marker gene panels
    df = pd.read_csv(args.multigene_table,sep='\t',index_col=0)
    df = df.loc[:,args.genes]
    gi = list(set(df.values.flatten().tolist()))

    # species that had hits
    gn = [v for v in refmap.values()]
    # flatten arrays
    gn = list(itertools.chain.from_iterable(gn))
    gn = list(set(gn))

    mgseq = {}
    for record in seq_util.iter_seq(args.pangenome_marker):
        nm = record[0][1:].split(' ')
        strand = nm[-1]
        nm = nm[0]
        mgseq[nm]= record[1].upper()
                
    # write          
    for i in gi:
        if i in refmap:
            ofn = os.path.join(args.output_dir, i+'.ref.fna')
            seq = {k: mgseq[k] for k in refmap[i]}
            tree_util.dict2fasta(seq, ofn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--filtered_blast', type=str, help='Filtered blast result path')
    parser.add_argument('--multigene_table', type=str, help='Multi gene table path')
    parser.add_argument('--pangenome_marker', type=str, help='Sequence file path to the pangenome marker genes')
    parser.add_argument('--genes', metavar='G', type=str, nargs='+', help='List of marker genes')
    parser.add_argument('--output_dir', type=str, help='Output file format')
    args = parser.parse_args()
    main(args)
