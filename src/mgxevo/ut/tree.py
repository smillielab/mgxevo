import sys
from . import seq as seq_util
import re, pickle,glob,os
import numpy as np
import pandas as pd
import subprocess

import random, string
def fread(file_path):
    import csv
    import pandas as pd
    sniffer = csv.Sniffer()
    data = open(file_path, "r").read(4096)
    delimiter = sniffer.sniff(data).delimiter
    return pd.read_csv(file_path, sep=delimiter)
# blast query to hits
def blast2dict(fn):
    df = pd.read_csv(fn ,sep='\t',header=None).set_axis('qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(), axis=1)
    x = df.qseqid.tolist()
    y = df.sseqid.tolist()
    d = {i:[] for i in x}
    for i in range(0,len(x)):
        d[x[i]].append(y[i])
    return(d)

def blast_pd2dict(df):
    x = df.iloc[:,0].tolist()
    y = df.iloc[:,1].tolist()
    d = {i:[] for i in x}
    for i in range(0,len(x)):
        d[x[i]].append(y[i])
    return(d)

# read BLAST files with nice headers
def read_blast(fn):
    return(pd.read_csv(fn ,sep='\t',header=None).set_axis('qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(), axis=1))

# parse qacct stats
def parse_qacct(job_script):
    x = subprocess.run(['qacct', '-j',job_script], stdout=subprocess.PIPE).stdout.decode('utf-8')
    x = x.split('==============================================================')[1:]
    x = {int(re.search(r'jobnumber\s+(\d+)\s+\n',xi).group(1)): {line.split()[0]:' '.join(line.split()[1:]) for line in xi.split('\n')[1:-1]} for xi in x}
    req_mem = {k: int(re.search('h_vmem=(\d+)g',xi['submit_cmd']).group(1))*1e9 for k, xi in x.items()}
    max_mem = {k: xi['maxvmem'] for k, xi in x.items()}

    scale = {}
    for k,v in max_mem.items():
        if re.search('M', v):
            scale[k] = 1e6*float(v[:-1])
            continue
        if re.search('G', v):
            scale[k] = 1e9*float(v[:-1])
            continue
        if re.search('K', v):
            scale[k] = 1e3*float(v[:-1])
            continue
        scale[k] = 1*float(v[:-1])
    num_jobs = len(x)
    exc_timelim = [str(k) for k,v in x.items() if re.search('h_rt limit', v['failed'])]
    exc_memlim = [str(k) for k,v in x.items() if scale[k] > req_mem[k]]
    output = ['From %s jobs:' %(num_jobs)]
    if(len(exc_timelim)>0):
        output += ['The following %s job(s) exceeded their time limits: %s' %(str(len(exc_timelim)), ', '.join(exc_timelim))]
    else:
        output += ['No job(s) exceeded their time limits']
    if(len(exc_memlim)>0):
        output += ['The following %s job(s) exceeded their memory limits: %s' %(str(len(exc_memlim)), ', '.join(exc_memlim))]
    else:
        output += ['No job(s) exceeded their memory limits']
    output = '\n\t'.join(output)
    print(output)
    cmds = readLines(x[list(x.keys())[0]]['cwd']+'/error/'+x[list(x.keys())[0]]['jobname'])[5:-3]
    cmds = gsub('^array\[\d+\]="|"$','',cmds)
    task2jobnum = {x[k]['taskid']: k for k in x.keys()}
    cmds = {task2jobnum[str(si)]: cmds[si-1] for si in range(1,len(cmds)+1)}
    return((exc_timelim, exc_memlim, x, cmds))

## 24 JULY
# instead of temp.fna use random string to prevent file clashes
def uniq(x):
    return list(set(x))

def gsub(pattern, replacement, string):
    if isinstance(string,str):
        return re.sub(pattern, replacement, string)
    elif isinstance(string,list):
        return [re.sub(pattern, replacement, x) for x in string]

def cbind(data):
    return pd.concat(data, axis=1)

def rbind(data):
    return pd.concat(data, axis=0)

def writeLines(x, file, append = False):
    if isinstance(x,str):
        x = [x]
    if append:
        with open(file,'a') as f:
            print('\n'.join(x),file=f)
    else:
        with open(file,'w') as f:
            print('\n'.join(x),file=f)

def readLines(ifn):
    with open(ifn,'r') as f:
        res = f.read().splitlines()
    
    if(len(res)==1):
        res = res[0]
    
    return res 

def table(data):
    if(type(data)=='list'):
        return pd.Series(data).value_counts()
    else:
        return pd.DataFrame(data).value_counts()

def getfiles(path, ext):
    ret = glob.glob(path+'*'+ext)
    ret = gsub(path,'',ret)
    ret = gsub(ext,'',ret)
    return ret

def unpickle(filepath, opt = 'rb'):
    ret = pickle.load(open(filepath,opt))
    return ret
    
def gsearch(pattern, string, invert=False):
    matches = string
    if not invert:
        for p in tuple(pattern):
            matches = [xq for xq in matches if re.search(p,xq) is not None]

    else:
        for p in tuple(pattern):
            matches = [xq for xq in matches if re.search(p,xq) is None]

    return matches

def reverse_complement(x):
    return x[::-1].translate(x.maketrans('acgtnACGTN-', 'tgcanTGCAN-'))

def charswap(original, locations, replacements):
    seq = np.array(list(original))
    seq[locations] = replacements
    seq = ''.join(seq)
    return seq
    
def fasta2np(fn):
    acgt = {'A':[1,0,0,0,0], 'C':[0,1,0,0,0], 'G':[0,0,1,0,0], 'T':[0,0,0,1,0], 'N':[0,0,0,0,1]}
    x = []
    i = []
    for record in seq_util.iter_seq(fn):
        i.append(record[0][1:])
        x.append([acgt[nt] if nt in acgt else [0,0,0,0,0] for nt in record[1].upper()])
    x = np.array(x)
    i = np.array(i)
    j = np.array(range(x.shape[1]))
    return {'x':x, 'i':i, 'j':j} 


def get_cons_genotypes(x, ind=None, seq=None):
    # return [1 x L] nucleotide sequence
    # x = alignment object (keys=x,i,j)
    # seq = sequence object (keys=fst, beg)
    acgt = np.array('A C G T N'.split())
    if ind is None:
        ind = range(x['x'].shape[0])
    if type(ind) != 'list':
        ind = list(ind)
    if seq is None:
        g1 = ''.join(acgt[list(x['x'][ind,:,:].sum(axis=0).argmax(axis=1))])
    else:
        u = x
        g1 = np.array([np.array([nt for nt in seq['fst']])]*u['x'].shape[0])
        # if empty mark with N
        nts = u['x'].argmax(axis=2)
        nts[u['x'].max(axis=2)==0] = 4
        g2 = acgt[nts]
        i = x['j'] - seq['beg']
        g1[:,i] = g2
        g1 = np.array([''.join(g) for g in g1])
        g1 = g1[ind]
    return g1

def get_genotypes(x, ind=None, seq=None):
    # return [1 x L] nucleotide sequence
    # x = alignment object (keys=x,i,j)
    # POSITION STARTS AT 1
    # seq = sequence object (keys=fst, beg)
    acgt = np.array('A C G T N'.split())
    if ind is None:
        ind = range(x['x'].shape[0])
    if type(ind) != 'list':
        ind = list(ind)
    if seq is None:
        g1 = ''.join(acgt[list(x['x'][ind,:,:].sum(axis=0).argmax(axis=1))])
    else:
        g1 = np.array([nt for nt in seq['fst']])
        g2 = acgt[list(x['x'][ind,:,:].sum(axis=0).argmax(axis=1))]
        i = x['j'] - seq['beg']
        g1[i] = g2
        g1 = ''.join(g1)
    
    return g1
    
def np2fasta(x, lab=None, seq=None, out='test.fst'):
    if lab is None:
        lab = x['i']
    fh = open(out, 'w')
    for i in range(x['x'].shape[0]):
        fh.write('>%s\n%s\n' %(lab[i], get_genotypes(x, ind=[i], seq=seq)))
    fh.close()

def np2dict(x, lab=None, seq=None):
    d = {}
    if lab is None:
        lab = x['i']
    # fh = open(out, 'w')
    for i in range(x['x'].shape[0]):
        d[lab[i]] = get_genotypes(x, ind=[i], seq=seq)
    return d

def np2dict2(x, lab=None, seq=None):
    d = {}
    if lab is None:
        lab = x['i']
    # fh = open(out, 'w')
    for i in range(x['x'].shape[0]):
        d[lab[i]] = get_cons_genotypes(x, ind=[i], seq=seq)
    return d

def org2seq(org):
    # read fasta sequences in a dictionary
    seq = {}
    for line in open('/broad/smillie-data/proj/mgxevo/map_reads/gtf/%s.gene_file.txt' %(org)):
        line = line.rstrip().split('\t')
        try:
            [contig, gene, beg, end, strand, fst] = line
        except:
            continue
        seq[gene] = {}
        seq[gene]['beg'] = int(beg)
        seq[gene]['end'] = int(end)
        seq[gene]['fst'] = fst
        seq[gene]['strand'] = strand
    
    return seq

def filtered_np_index(species, gi):
    samples = readLines('/broad/smillie-data/proj/mgxevo/samples.txt')
    u = unpickle('/broad/smillie-data/proj/mgxevo/gwas/%s/%s.majority_snp.pickle' %(species, gi))
    v = unpickle('/broad/smillie-data/proj/mgxevo/map_reads/temp/%s/%s.aln.pickle' %(species,gi))
    # convert names to sample indicies
    iu = [samples.index(fix_suffix(k)) for k in u.keys()]
    ix = list(set(iu).intersection(v['i']))
    # iterate to see which indicies in v['i'] we use
    use  = [si in ix for si in v['i']]
    return use

def filtered_np_index2(aln, msnp, samples='/broad/smillie-data/proj/mgxevo/samples.txt'):
    # filtered_np_index without hardcoded paths (csmillie)
    # input arguments:
    # aln = alignment (pickle)
    # msnp = majority snp (pickle)
    # samples = list of samples
    # read input data
    samples = readLines(samples)
    u = unpickle(msnp)
    v = unpickle(aln)
    # convert names to sample indicies
    iu = [samples.index(fix_suffix(k)) for k in u.keys()]
    ix = list(set(iu).intersection(v['i']))
    # iterate to see which indicies in v['i'] we use
    use  = [si in ix for si in v['i']]
    return use

def dict2fasta(d,ofn):
    with open(ofn, 'w') as f: [f.write(">" + i + "\n" + d[i] + "\n") for i in d.keys()]

def prepend(string, ifn):
    with open(ifn, 'r') as original: data = original.read()
    with open(ifn, 'w') as modified: modified.write(string + "\n" + data)

def dict2np(d):
    tmp_file = id_generator(64)+'.fna'
    dict2fasta(d,tmp_file)
    ret = fasta2np(tmp_file)
    # remove tempfile
    os.system('rm %s*' %(tmp_file))
    return ret
    
def add_fasta2np(x, lab, seq, fna):
    tfn1 = id_generator(64)+'.fna'
    tfn2 = id_generator(64)+'.fna'
    tfn3 = id_generator(64)+'.fna'

    np2fasta(x, lab=lab, seq=seq, out=tfn1)
    os.system('cat %s %s > %s' %(tfn1, fna, tfn2))
    os.system('~/bin/mafft %s > %s' %(tfn2, tfn3))
    res = fasta2np(tfn3)
    j = x['j'] - seq['beg']
    res['x'] = res['x'][:,j,:]
    res['j'] = x['j']
    
    # remove temp files
    os.system('rm %s*' %(tfn1))
    os.system('rm %s*' %(tfn2))
    os.system('rm %s*' %(tfn3))
    return res

def add_fasta2tree(tree_fas, fna, ofn):
    # input: XXX.tree (assumes XXX.tree.fas is in same dir)
    tfn1 = id_generator()
    tfn2 = id_generator()

    os.system('cat %s %s > %s.fna' %(tree_fas, fna, tfn1))
    os.system('~/bin/mafft %s.fna > %s.mafft' %(tfn1,tfn2))
    # new alignments are stored in temp.mafft
    prefix = '____iqtree.temp.' + id_generator()
    os.system('python ~/bin/run_gblocks.py --aln %s.mafft --lax' %(tfn2))
    os.system('~/bin/iqtree2 -s %s.mafft.gb -m GTR+G4 -alrt 1000 --prefix %s -fast --seed 499244' %(tfn2,prefix))

    # move tree file
    os.system('mv %s.treefile %s' %(prefix, ofn))
    
    # copy unique sequences
    if os.path.isfile(prefix+'.uniqueseq.phy'):
        os.system('mv %s.uniqueseq.phy %s.uniq' %(prefix, ofn))

    # move log
    os.system('mv %s.log %s.log' %(prefix, ofn))
        
    # remove extra IQ tree files
    if len(glob.glob('%s*' %(prefix)))>0:
        os.system('rm %s*' %(prefix))
    
    # remove temp files
    os.system('rm %s*' %(tfn1))
    os.system('rm %s*' %(tfn2))
    
def fix_suffix(samples):
    return gsub('_(UC|CD|HC|IBDU|HR)$','',samples)

def id_generator(N=64):
    return ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in range(N))

def match(a,b):
    return [ b.index(x) if x in b else None for x in a ]
