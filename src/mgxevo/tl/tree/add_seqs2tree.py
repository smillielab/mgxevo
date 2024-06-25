import sys
import argparse
import re
import pickle
import glob
import os
import shutil
import logging
import subprocess
from pathlib import Path
from Bio import Phylo

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

import ut
from ut import seq as seq_util
from ut import tree as tree_util

def get_reference_sequence(gtf_file):
    """
    Get reference sequence from a gtf file.
    """
    fst = {}
    for line in open(gtf_file):
        contig, gene, beg, end, strand, seq = line.rstrip().split('\t')
        fst[gene] = {'beg':int(beg), 'end':int(end),'fst':seq, 'strand':strand}
    return fst

def write_to_file(file_name, content):
    with open(file_name, 'w') as f:
        f.write(content)

def get_new_sequences(seq_file):
    name = {}
    new_seqs = {}
    for record in seq_util.iter_seq(seq_file):
        nm = record[0][1:]
        new_nm = re.sub('[^a-zA-Z0-9]','',nm)
        name[new_nm] = nm
        new_seqs[new_nm]= record[1].upper()
    return name, new_seqs

def main(args):
    
    # get reference sequence
    fst = get_reference_sequence(args.gtf)

    seq = f'>REF\n{fst[args.name]["fst"]}'

    prefix = f'{args.tmp_dir}/__usher.{args.name}'
    prefix2 = f'{prefix}/__usher.{args.name}'

    os.makedirs(prefix, exist_ok=True)

    write_to_file(f'{prefix2}.1.fna', '\n'.join([seq,open(args.fasta,'r').read()]))

    ut.run_cmd(f'mafft --auto {prefix2}.1.fna > {prefix2}.1.aln.fna', silent=True)
    ut.run_cmd(f'faToVcf -ref=REF {prefix2}.1.aln.fna {prefix2}.1.vcf', silent=True)

    # the new sequences can't contain any characters so strip them but keep track
    name, new_seqs = get_new_sequences(args.seq)

    write_to_file(f'{prefix2}.ref.fna', seq)

    ofn = f'{prefix2}.2.fna'
    tree_util.dict2fasta(new_seqs,f'{ofn}temp')

    # align
    logging.info('Aligning new sequences to tree.')
    cmd =f'mafft --auto --keeplength --addfragments {ofn}temp {prefix2}.ref.fna > {ofn}'
    ut.run_cmd(cmd, silent=True)

    ut.run_cmd(f'faToVcf -ref=REF {ofn} {prefix2}.2.vcf', silent=True)
    logging.info('Initialising tree.')
    ut.run_cmd(f'usher -t {args.tree} -v {prefix2}.1.vcf -o {prefix2}.pb -d {prefix}', silent=True)
    logging.info('Appending new sequences to tree.')

    if args.save:
        ut.run_cmd(f'usher -v {prefix2}.2.vcf -i {prefix2}.pb -u -d {prefix} -o {prefix2}.updated.pb', silent=True)
    else:
        ut.run_cmd(f'usher -v {prefix2}.2.vcf -i {prefix2}.pb -u -d {prefix}', silent=True)

    # fix any branches that haven't been assigned a length (occurs when the tip is assigned at a node)
    ifn = f'{prefix}/uncondensed-final-tree.nh'
    br = Phylo.read(ifn,'newick').depths()

    bl=[f.branch_length for f in br.keys() if f.branch_length is not None]

    r=[f.name for f in br.keys() if f.branch_length is None]

    tree = tree_util.readLines(ifn)
    if r:
        missing_samples = "\n".join(r)
        logging.warning(f'Following samples are missing branch lengths (defaulting to 1e-8): {missing_samples}')
        for rq in r:
            tree = tree_util.gsub(rq,f'{rq}:0.00000001',tree)

    tree_util.writeLines(tree,args.output)

    # tidy up
    if args.save:
        
        # shutil.move(f'{prefix2}.vcf', f'{args.output}.vcf')
        shutil.move(f'{prefix2}.pb', f'{args.output}.pb')
        shutil.move(f'{prefix2}.updated.pb', f'{args.output}.updated.pb')

    shutil.rmtree(prefix, ignore_errors=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', help='tree (newick)')
    parser.add_argument('--name', help='reference name')
    parser.add_argument('--fasta', help='tree (fasta)')
    parser.add_argument('--gtf',  help='gtf file that contains ref seq for name')
    parser.add_argument('--seq', help='seqs to be added to tree (fasta)')
    parser.add_argument('--output', help='output tree (newick)')
    parser.add_argument('--tmp_dir', help='tmp dir', default='./tmp')
    parser.add_argument('--save', help='save USHER pb file?', default=True, action='store_true')

    args = parser.parse_args()
    main(args)
    