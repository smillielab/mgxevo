import re, pickle, glob, os
import sys
import argparse
from pathlib import Path
import shutil
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))

import ut
from ut import seq as seq_util
from ut import tree as tree_util

parser = argparse.ArgumentParser()
parser.add_argument('--name', help='Gene name')
parser.add_argument('--tree_dir', help='Directory containing the tree')
parser.add_argument('--output_dir', help='Directory to write the output')
parser.add_argument('--multi', help='Whether to use multi or single', action='store_true')
parser.add_argument('--refseq_path', help='Path to the reference sequence', required=True)
parser.add_argument('--multigene_table', help='Path to the multigene table, required if running multigene analysis', required=False)
parser.add_argument('--genes', metavar='G', type=str, nargs='+', help='List of marker genes to be processed, required if running multigene analysis', required=False)
args = parser.parse_args()

refs = {}
if args.refseq_path is None or not os.path.exists(args.refseq_path):
    raise Exception(f'`args.refseq_path` does not exist or not provided')
for records in seq_util.iter_seq(args.refseq_path):
    nm = records[0][1:]
    refs[nm] = records[1].upper()

if args.multi:
    if args.multigene_table is None or not os.path.exists(args.multigene_table):
        raise Exception(f'`args.multigene_table` does not exist or not provided')
    if args.genes is None:
        raise Exception(f'`args.genes` not provided')
    multigene = pd.read_table(args.multigene_table, index_col=0)
    for gn in multigene.index:
        refs[gn] = ''.join([refs[multigene.loc[gn][v]] for v in args.genes])
    tmp = tree_util.id_generator(8)
else:
    tmp = args.name

# make database
path = f'{args.tree_dir}/{args.name}'
os.makedirs(path, exist_ok=True)

ofn = f'{path}/{tmp}.ref.fa' # the ref seq for the tree
tree_util.writeLines(['>REF',refs[args.name]], ofn)

cmd = f'cd {path}; samtools faidx {tmp}.ref.fa'
ut.run_cmd(cmd, silent=True)

# PB to VCF
ifn = f'{args.tree_dir}/{args.name}.ref.tree.updated.pb'
cmd = f'matUtils extract -i {ifn} -d {path} -v {tmp}.vcf'
ut.run_cmd(cmd, silent=True)

# VCF to fasta
ifn1 = os.path.realpath(f'{path}/{tmp}.ref.fa') # the ref seq for the tree
ifn2 = os.path.realpath(f'{path}/{tmp}.vcf') # the VCF for the tree

cmd = f'cd {path}; vcf2fasta -f {ifn1} {ifn2}'
ut.run_cmd(cmd, silent=True)

# # tidy up
# ifn1 = f'{path}/{tmp}.ref.fa' # the ref seq for the tree
# ifn2 = f'{path}/{tmp}.vcf' # the VCF for the tree
# ifn3 = f'{path}/{tmp}.ref.fa.fai' # index file
# cmd = f'rm {ifn1} {ifn2} {ifn3}'
# ut.run_cmd(cmd, silent=True)

# ofn = f'{args.name}.usher.fst'
# cmd = f'cd {path}; cat *.fa > {ofn}; rm *.fa'
# ut.run_cmd(cmd, silent=True)

# # fix seq names:
# ifn = f'{path}/{args.name}.usher.fst'
# ofn = f'{args.output_dir}/{args.name}.usher.fst'
# seq = seq_util.read_fst(ifn)
# seq  = {re.sub('\\_REF:\\d+$','',k): v for k,v in seq.items()}
# tree_util.dict2fasta(seq, ofn)

# shutil.rmtree(path, ignore_errors=True)