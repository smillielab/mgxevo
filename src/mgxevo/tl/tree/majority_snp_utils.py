import sys, os, pickle, re, glob
import pandas as pd
import numpy as np 
import scipy.stats as stats
from Bio import Phylo
from Bio import SeqIO
import subprocess
import logging
import shutil
from pathlib import Path
import argparse
from collections import defaultdict

from .. import tree
from ... import ut
from ...ut import tree as tree_util
from ...ut import seq as seq_util
from ...data import metadata_path

def fix_tree(ofn, thres=75):
    logging.debug(f"Fixing tree {ofn}")
    ut.run_cmd(f'Rscript {tree.script_path}/tools/fix_tree.r {ofn}.treefile {thres} {ofn}.treefile', silent=True)
    

# [FASTA files] + output loc --> multigene phylogeny 
def multigene_phylo(fasta_path, gene_ids, output, fasta_ext = 'fas.processed.gb', fast=True, tmp_dir='./tmp', name=None, seed=None, threads=None):
    # check if exist
    fasta_files = [Path(fasta_path) / f"{gene_id}.{fasta_ext}" for gene_id in gene_ids]
    if not all(fn.exists() and fn.stat().st_size > 0 for fn in fasta_files):
        bad = [fn for fn in fasta_files if not fn.exists() or fn.stat().st_size == 0]
        missing_fst = '\n'.join(map(str, bad))
        raise NameError(f"The following FASTAs do not exist: {missing_fst}")

    # add random id to tmp_dir to avoid conflicts between parallel runs
    tmp_dir = Path(tmp_dir) / tree_util.id_generator(8)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # nexus file
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if name is None:    
        name = gene_ids[0].rsplit('_', 1)[0]
        gi = gene_ids[0].rsplit('_', 1)[0]
    else:
        gi = name
        
    if threads is None:
        threads = 'AUTO'

    ofn = output_path / f'{gi}.nex'
    nex = '#nexus\nbegin sets;\n' + '\n'.join([f'    charset part{i+1} = {fn}: *;' for i, fn in enumerate(fasta_files)]) +  '\nend;\n'
    ofn.write_text(nex)

    ofn1 = Path(tmp_dir) / f'___iqtree___{gi}'
    cmd = ['iqtree2', '-p', str(ofn), '-m', 'GTR+G4' if fast else 'MFP', '-alrt', '1000', '--prefix', str(ofn1), '-nt', str(threads)]
    if seed is not None:
        cmd += ['-seed ', str(seed)]

    if fast:
        cmd.append('-fast')
    ut.run_cmd(cmd, silent=True)

    # fix tree
    fix_tree(ofn1)
    
    # move tree file
    shutil.move(str(ofn1) + '.treefile', str(output_path))
    (output_path / f'{ofn1.name}.treefile').rename(output_path / f'{name}.tree')

    # copy unique sequences
    uniq_file = ofn1.with_suffix('.uniqueseq.phy')
    if uniq_file.exists():
        shutil.move(str(uniq_file), output_path.with_suffix(f'.{gi}.uniq'))

    # move log
    shutil.move(str(ofn1) + '.log', output_path.with_suffix(f'.{gi}.log'))

    # remove other files
    for file in ofn1.parent.glob(f'{ofn1.name}*'):
        file.unlink()
    
# Get a phylogeny for a given species and mean depth filtering threshold
def get_phylogeny_from_pickle(pickle_file, gene_file, gi=None, dx_file=None, output_dir='./write/', fast=True, seed=None, temp_dir='./tmp', threads=None):
    
    def run_gblocks(ofn):
        cmd = ['python', f'{tree.script_path}/tools/run_gblocks.py', '--aln', str(ofn), '--lax']
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            logging.error(f"Gblocks failed with error: {result.stderr.decode()}")
            raise RuntimeError("Gblocks execution failed.")

    def run_iqtree(ifn1, ofn1, fast, seed, threads):
        cmd = ['iqtree2', '-s', str(ifn1), '-m', 'GTR+G4' if fast else 'MFP', '-alrt', '1000', '--prefix', str(ofn1), '-nt', str(threads)]
        if fast:
            cmd.append('-fast')
        if seed is not None:
            cmd.extend(['-seed', str(seed)])
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            logging.error(f"IQTREE failed with error: {result.stderr.decode()}")
            raise RuntimeError("IQTREE execution failed.")

    def move_files(ofn1, output_path):
        # move tree file
        ofn1.with_suffix('.treefile').rename(output_path.with_suffix('.tree'))
        # copy unique sequences
        uniq_file = ofn1.with_suffix('.uniqueseq.phy')
        if uniq_file.exists():
            uniq_file.rename(output_path.with_suffix('.uniq'))
        # move log
        log_file = ofn1.with_suffix('.log')
        if log_file.exists():
            log_file.rename(output_path.with_suffix('.log'))
        # remove other files
        for file in Path(ofn1.parent).glob(f'{ofn1.name}*'):
            file.unlink()
    
    if threads is None:
        threads = 'AUTO'
        
    temp_dir = Path(temp_dir) / tree_util.id_generator(8)
    
    if gi is None:
        gi = os.path.basename(pickle_file).replace('.aln.pickle','')
    if dx_file is None:
        dx_file = metadata_path
    
    ret = majority_snps(pickle_file, gene_file, gi, dx_file=dx_file, min_depth=3, threads=threads, tmp_dir=temp_dir)
    output_path = Path(output_dir) / f'{gi}'
    ofn = output_path.with_suffix('.fas')
    tree_util.dict2fasta(ret, str(ofn))
    run_gblocks(ofn)
    ifn1 = ofn.with_suffix('.fas.gb')
    ofn1 = Path(temp_dir) / f'___iqtree___{gi}'
    run_iqtree(ifn1, ofn1, fast, seed, threads)
    fix_tree(ofn1)
    move_files(ofn1, output_path)
    

# Return a dictionary of the majority SNPs for each sample (picking the subject)
def majority_snps(pickle_file: str, gene_file: str , gi=None, dx_file=None, min_depth=3, max_std_err=3, phylo_filter=True, tmp_dir='./tmp', threads=None):
    """
    Returns a dictionary of the majority SNPs for each sample (picking the subject).

    Args:
        pickle_file (str): Path to pickle file of alignment.
        gene_file (str): Gene file in format of table (.txt).
        gi: Gut genome id.
        dx_file (str): Sample table with diagnosis.
        min_depth (int, optional): Minimum depth. Defaults to 3.
        max_std_err (int, optional): Maximum standard error. Defaults to 3.
        phylo_filter (bool, optional): Whether to apply phylogenetic filter. Defaults to True.

    Returns:
        dict: Dictionary of the majority SNPs for each sample.
    """
    
    def run_iqtree(input_file, prefix, threads, seed=None):
        IQTREE_CMD = 'iqtree2 -s {input} -m GTR+G4 -alrt 1000 --prefix {prefix} -fast -quiet -seed 12345 -nt {threads}'
        cmd = IQTREE_CMD.format(input=input_file, prefix=prefix, threads=threads)
        result = subprocess.run(cmd, shell=True, capture_output=True)
        if result.returncode != 0:
            logging.error(f"IQTREE failed with error: {result.stderr.decode()}")
            raise RuntimeError("IQTREE execution failed.")
        # Move tree file
        prefix.with_suffix('.treefile').rename(input_file.with_suffix('.fas.tree'))
        return result

    def cleanup_files(tmp_dir, gi):
        ofn = tmp_dir / f'{gi}.fas'
        files_to_remove = list(tmp_dir.glob(f'___iqtree___*')) + [ofn]
        for file in files_to_remove:
            if file.exists():
                file.unlink()

    def process_tree(ofn):
        BOOTSTRAP_THRESHOLD = 75
        tree = Phylo.read(f'{ofn}.tree', 'newick')
        tree.collapse_all(lambda c: c.confidence is not None and c.confidence < BOOTSTRAP_THRESHOLD)
        tree.root_at_midpoint()
        return tree

    def filter_samples(tree, max_std_err, indsnp):
        depth1 = tree.depths()
        depths = {str(key):depth1[key] for key in depth1.keys() if str(key) != 'Clade'}
        keys, vals = zip(*depths.items())
        z = abs(stats.zscore(vals))<max_std_err
        if all(np.isnan(stats.zscore(vals))):
            z = np.logical_not(z)
        newmap = dict(zip(keys,z))
        [indsnp.pop(k, None) for k,v in newmap.items() if not v]

    # Check if files exist
    if not os.path.isfile(pickle_file):
        raise FileNotFoundError(f"Pickle file {pickle_file} does not exist.")
    if not os.path.isfile(gene_file):
        raise FileNotFoundError(f"Gene file {gene_file} does not exist.")
    if not os.path.isfile(dx_file):
        raise FileNotFoundError(f"Sample table file {dx_file} does not exist.")
    if gi is None:
        gi = os.path.basename(pickle_file).replace('.aln.pickle','')
    
    # Process gene file
    fst = {}
    with open(gene_file) as file:
        for line in file:
            contig, gene, beg, end, strand, seq = line.rstrip().split('\t')
            fst[gene] = {'beg':int(beg), 'end':int(end),'fst':seq, 'strand':strand}

    # Load alignments of sample to a gut_genome
    with open(pickle_file, 'rb') as f:
        u = pickle.load(f)

    # Check if gene file is empty
    if 0 in u['x'].shape:
        return None

    # Filter samples with mean depth < threshold
    i =  u['x'].sum(axis=2).mean(axis=1) > min_depth
    u['x'] = u['x'][i,:,:]
    u['i'] = u['i'][i]
    
    # Load diagnosis data
    dx = pd.read_csv(dx_file, sep='\t').set_index('Sample ID')
    dx = dx.iloc[u['i']]
    dx['depth'] = u['x'].sum(axis=2).mean(axis=1)

    # Get index of maximum depth for each subject
    ix = dx.loc[dx.groupby("Subject")["depth"].idxmax()].index.tolist()
    ix = dx.index.get_indexer_for(ix)
    ix = ix[np.argsort(u['i'][ix])]

    u['x'] = u['x'][ix,:,:]
    u['i'] = u['i'][ix]

    # Create names for each sample
    names = dx.index + '_' + dx['Diagnosis (UC, CD, IBDU, HC)']
    names = names.iloc[ix]

    # Get majority SNP and mark entries where there's no max with N
    nts = u['x'].argmax(axis=2)
    nts[u['x'].max(axis=2)==0] = 4

    # Convert numpy object to genotypes
    seq = fst[gi]
    acgt = np.array(['A', 'C', 'G', 'T', 'N'])
    ind = range(u['x'].shape[0])
    j = u['j']-seq['beg']

    g = np.array([np.array([nt for nt in seq['fst']])]*u['x'].shape[0])
    g[:,j] = acgt[nts]

    # Note: output will have + orientation
    if seq['strand']=='-':
        X  = {names.iloc[i]: tree_util.reverse_complement(''.join(g[i,:])) for i in range(g.shape[0])}
    else:
        X  = {names.iloc[i]: ''.join(g[i,:]) for i in range(g.shape[0])}

    # Filter samples based on diagnosis
    df = dx.iloc[ix,].reset_index().set_index('Subject').loc[dx.Subject.unique(),]
    indsnp = {k:X[k] for k in df['Sample ID']+'_'+df['Diagnosis (UC, CD, IBDU, HC)']}

    if len(indsnp)<4:
        return indsnp
    
    if phylo_filter:
        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        ofn = tmp_dir / f'{gi}.fas'
        tree_util.dict2fasta(indsnp, ofn)
        run_iqtree(ofn, tmp_dir / f'___iqtree___{gi}', threads)
        tree = process_tree(ofn)
        filter_samples(tree, max_std_err, indsnp)
        cleanup_files(tmp_dir, gi)

    return indsnp


def process_fna_gb(fasta_path, gene_ids, ext = 'fas.gb', name=None, output_dir = None):
    u = {gid: seq_util.read_fst(f'{fasta_path}/{gid}.{ext}') for gid in gene_ids}

    pt = [list({len(v) for v in u[gi].values()})[0] for gi in gene_ids]
    pdn = {gi: list({len(v) for v in u[gi].values()})[0] for gi in gene_ids}


    beg = [1]+list(np.cumsum(pt)[:-1]+1)
    end = list(np.cumsum(pt))
    nex = '#nexus\nbegin sets;\n' + '\n'.join(['    charset part%s = %s-%s;' %(i+1,beg[i],end[i]) for i,fn in enumerate(pt)]) +  '\nend;\n'

    nms = []

    nms_col = sum([list(u[gi].keys()) for gi in gene_ids],[])
    for i in range(len(nms_col)):
        val = nms_col[i]
        if val not in nms:
                nms += [val]

    concat = {sample: re.sub('-','N',''.join([u[gi].get(sample, '-'*pdn[gi]) for gi in gene_ids])) for sample in nms}
    
    if output_dir is None:
        output_dir = fasta_path
    
    # output file paths
    ofn_fas = f'{output_dir}/{name}.conaln.fas'
    ofn_nex = f'{output_dir}/{name}.conaln.nex'
    
    # write FASTA
    out = ''.join(['>%s\n%s\n' %(nm, concat[nm]) for nm in nms])
    with open(ofn_fas,'w') as f:
        print(out,file=f)
        
    # write NEXUS file
    with open(ofn_nex,'w') as f:
        print(nex,file=f)