import os
import sys
import copy
import pickle
import argparse
import numpy as np
import pandas as pd
import scipy.stats
import sklearn
import sklearn.cluster
import ete3
import re

from pathlib import Path
sys.path.append(str(Path(__file__).resolve()))
from tools.StrainFinder import *

def disease_index():
    return dx.loc[np.array(samples)[x['gene']['i']], 'Diagnosis (UC, CD, IBDU, HC)']


# input arguments
# ---------------
parser = argparse.ArgumentParser()
parser.add_argument('--genome', help='Genome id')
parser.add_argument('--aln', help='Multigene alignment file')
parser.add_argument('--tree', help='Multigene tree file')
parser.add_argument('--samples', help='Samples file', default='samples.txt')
parser.add_argument('--metadata', help='Metadata file')
parser.add_argument('--final_tips', help='Final tips file', default='strainfinder.final_tips.txt')
parser.add_argument('--output', help='output file (.pickle)')
parser.add_argument('--use_phylo', help='only use samples that appear in the phylogeny', default = '1')
args = parser.parse_args()

# load data
# ---------

# read samples
samples = [line.rstrip() for line in open(args.samples)]

# read iHMP metadata
dx = pd.read_table(args.metadata, header=0, index_col=0)

# load alignment
x = {}
if os.path.exists(args.aln):
    x['gene'] = pickle.load(open(args.aln, 'rb'))
    print('%.2f%% N bases' %(100.*x['gene']['x'][:,:,4].sum() / x['gene']['x'].sum()))
    x['gene']['x'] = x['gene']['x'][:,:,:4]
    print('%.2f total bases' %(x['gene']['x'].sum()))
else:
    quit(args.aln)

# get samples and disease index
si = np.array(samples)[x['gene']['i']]
di = disease_index()

# get samples to use
# ------------------
keep = [re.sub('_[^_]*$', '', li) for li in ete3.Tree(args.tree, format=1).get_leaf_names()]

# load consensus
# --------------

con = {}

def binvec(a):
    a = np.where(a == max(a))[0]
    a = random.choice(a)
    b = [0,0,0,0]
    b[a] = 1
    return np.array(b)

# read individual strains
for line in open(args.final_tips):
    if line.startswith('name'):
        continue
    line = line.rstrip().split()
    [name, node, coef, padj, sas] = line

    if name in args.aln:
        sas = [re.sub('_[^_]*$', '', sa) for sa in sas.split(';')]
        assert all([sa in keep for sa in sas])
        i = np.array([sa in sas for sa in si])
        print('estimating strain from %d samples' %(sum(i)))
        y = x['gene']['x'][i,:,:]
        y = np.apply_along_axis(binvec, 2, y)
        y = y.sum(axis=0, keepdims=True)
        y = np.apply_along_axis(binvec, 2, y)
        con[node] = y

# concatenate strains
nodes = sorted(con.keys())
y = np.concatenate([con[i] for i in nodes])

# strain finder
# -------------

# setup em object
print('reading input args', flush=True)
args2 = parse_args()
args2.msg = True

print('subsetting samples', flush=True)
assert(len(si) == x['gene']['x'].shape[0])
if bool(int(args.use_phylo)):
    i = np.array([sa in keep for sa in si])
else:
    i = np.array([sa in si for sa in si])

print('setting up object', flush=True)
x['gene']['x'] = x['gene']['x'][i,:,:]
si = si[i]
di = di[i]
data = Data(x=x['gene']['x'])
assert len(si) == x['gene']['x'].shape[0]

acgt = np.array(['a','c','g','t'])
seqs = [''.join(acgt[y[i,:,:].argmax(axis=1)]) for i in range(y.shape[0])]
out = open(f'{args.output}/%s.sf_cons.tab' %(args.genome), 'w')
for i,seq in enumerate(seqs):
    out.write('%s\t%s\n' %(nodes[i], seq))
out.close()

args2.data = data
em = load_em(args2)
nts = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])

n = y.shape[0] + 1

print('setting up estimate', flush=True)
estimate = Estimate(em.data, n=n, random=False, random_k=n, fix_p=y.shape[0])
estimate.p[:y.shape[0],:,:] = y[:,:,:4]

print('estimate strains', flush=True)
for it in range(4):
    p = copy.deepcopy(estimate.p)
    estimate.max_loglik_z()
    estimate.max_loglik_p()
    out = open('%s.%s' %(args.output, it), 'w')
    for j in range(estimate.z.shape[0]):
        out.write('%s_%s\t%s\n' %(si[j], di[j], '\t'.join(map(str, estimate.z[j,:]))))
    out.close()
    for i in range(estimate.p.shape[0]):
        print('iter %s; strain %s; s0 %s' %(it, i, (estimate.p[i,:,:] * y[:,:,:4]).sum(axis=2).mean()))
        print('iter %s; strain %s; si %s' %(it, i, (estimate.p[i,:,:] * p[i,:,:]).sum(axis=1).mean()))
