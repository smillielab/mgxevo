import sys
import argparse
import re
import pickle
import glob
import os
import shutil
import logging
import subprocess
import pandas as pd
from pathlib import Path
from Bio import Phylo

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

import ut
from ut import seq as seq_util
from ut import tree as tree_util

def get_reference_sequence_and_clean_names(args):
    tmp = tree_util.id_generator(8)
    tmp = f'{args.tmp_dir}/{tmp}'
    refs = {}
    if args.refseq_path is None or not os.path.exists(args.refseq_path):
        raise Exception(f'`args.refseq_path` does not exist or not provided')
    for records in seq_util.iter_seq(args.refseq_path):
        nm = records[0][1:]
        refs[nm] = records[1].upper()

    if args.multigene_table is None or not os.path.exists(args.multigene_table):
        raise Exception(f'`args.multigene_table` does not exist or not provided')
    if args.genes is None:
        raise Exception(f'`args.genes` not provided')
    multigene = pd.read_table(args.multigene_table, index_col=0)
    for gn in multigene.index:
        refs[gn] = ''.join([refs[multigene.loc[gn][v]] for v in args.genes])

    return refs, tmp

def get_sequence_lengths(fas):
    seqlen = []
    tree = {}
    for record in seq_util.iter_seq(fas):
        nm = record[0][1:]
        tree[nm]= record[1].upper()
        seqlen.append(len(record[1].upper()))
    return seqlen, tree

def write_ref_seq_to_file(refs, gn, tmp, seqlen):
    tree_util.writeLines(['>REF',refs[gn]], f'{tmp}.ref.fna')
    if (set([len(refs[gn])]) ^ set(seqlen) != set()):
        logging.warning(f'Length of reference sequence ({len(refs[gn])}) is not equal to length of sequences in the alignment ({seqlen}).')
        tree_util.writeLines(['>REF',refs[gn]] + tree_util.readLines(f'{tmp}.tree.fna'), f'{tmp}.tree.fna')
        ut.run_cmd(f'mafft --auto --keeplength {tmp}.tree.fna > {tmp}.tree.aln.fna', silent=True)
    else:
        tree_util.writeLines(['>REF',refs[gn]] + tree_util.readLines(f'{tmp}.tree.fna'), f'{tmp}.aln.fna')

def get_new_sequences(seq):
    name = {}
    new_seqs = {}
    for record in seq_util.iter_seq(seq):
        nm = record[0][1:]
        new_nm = re.sub('[^a-zA-Z0-9]','',nm)
        name[new_nm] = nm
        new_seqs[new_nm]= record[1].upper()
    return name, new_seqs

def fix_tree_that_have_no_lengths(ifn, out):
    nwk = tree_util.readLines(ifn)
    br = Phylo.read(ifn,'newick').depths()
    bl=[f.branch_length for f in br.keys() if f.branch_length is not None]
    r=[f.name for f in br.keys() if f.branch_length is None]
    if(len(r)):
        print('\nWARNING: Following samples are missing branch lengths (defaulting to 1e-8):\n'+'\n'.join(r)+'\n')
        for rq in r:
            nwk = tree_util.gsub(rq,f'{rq}:%s' %('0.00000001'),nwk)
    tree_util.writeLines(nwk, out)

def main(args):
    
    gn  = args.name
    seq = args.seq
    fas = args.fasta
    tre = args.tree
    out = args.output
    out_dir = os.path.dirname(out)

    refs, tmp = get_reference_sequence_and_clean_names(args)
    seqlen, tree = get_sequence_lengths(fas)
    tree_util.dict2fasta(tree, f'{tmp}.tree.fna')

    write_ref_seq_to_file(refs, gn, tmp, seqlen)

    ut.run_cmd(f'faToVcf -ref=REF {tmp}.aln.fna {tmp}.tree.vcf', silent=True)
    ut.run_cmd(f'usher -t {tre} -v {tmp}.tree.vcf -o {tmp}.pb -d {tmp}', silent=True)

    _, new_seqs = get_new_sequences(seq)

    tree_util.dict2fasta(new_seqs, f'{tmp}.newseqs.fna')

    ut.run_cmd(f'mafft --auto --keeplength --addfragments {tmp}.newseqs.fna {tmp}.ref.fna > {tmp}.newseqs.aln.fna', silent=True)
    ut.run_cmd(f'faToVcf -ref=REF {tmp}.newseqs.aln.fna {tmp}.newseqs.vcf', silent=True)
    ut.run_cmd(f'usher -v {tmp}.newseqs.vcf -i {tmp}.pb -u -o {tmp}.update.pb -d {tmp}/', silent=True)

    ifn = f'{tmp}/uncondensed-final-tree.nh'
    fix_tree_that_have_no_lengths(ifn, out)

    shutil.move(f'{tmp}.pb', f'{out_dir}/{gn}.pb')
    shutil.move(f'{tmp}.update.pb', f'{out_dir}/{gn}.ref.tree.updated.pb')
    shutil.rmtree(f'{tmp}')
    for f in glob.glob(f'{tmp}*'):
        os.remove(f)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', help='tree (newick)')
    parser.add_argument('--name', help='reference name')
    parser.add_argument('--fasta', help='tree (fasta)')
    parser.add_argument('--seq', help='seqs to be added to tree (fasta)')
    parser.add_argument('--refseq_path', help='Path to the reference sequence', required=True)
    parser.add_argument('--multigene_table', help='Path to the multigene table, required if running multigene analysis', required=False)
    parser.add_argument('--genes', metavar='G', type=str, nargs='+', help='List of marker genes to be processed', required=True)
    parser.add_argument('--tmp_dir', help='tmp dir', default='./tmp')
    parser.add_argument('--output', help='output tree (newick)')

    args = parser.parse_args()
    
    main(args)