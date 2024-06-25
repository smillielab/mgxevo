import numpy as np
import os, pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--strains", help="Path to genomes.use.txt file")
parser.add_argument("--aln_pickle_dir", help="Path to alignment pickle directory")
parser.add_argument("--multigene_table", help="Path to ids.multigene.txt file")
parser.add_argument("--genes", metavar='G', type=str, nargs='+',
                        help='a list of marker genes to be processed')
parser.add_argument("--output_dir", help="Path to output directory", default='.')
args = parser.parse_args()

with open(args.strains, 'r') as f:
    G = [line.rstrip() for line in f]

g = 0
with open(args.multigene_table, 'r') as f:
    header = f.readline().rstrip().split()
    for line in f:
        line = line.rstrip().split()
        genome = line[0]
        genes = [gene for gene, header_item in zip(line[1:], header[1:]) if header_item in args.genes]
        if genome not in G:
            continue

        if os.path.exists(os.path.join(args.output_dir, '%s.multigene.aln.pickle' %(genome))):
            print('skipping %s' %(genome))
            continue
        
        alignments = []
        for gene in genes:
            with open(os.path.join(args.aln_pickle_dir, '%s.aln.pickle' %(gene)), 'rb') as f:
                alignments.append(pickle.load(f))
        
        # merge alignments
        i = sorted(set(sum([list(aln['i']) for aln in alignments], [])))
        j = sum([list(aln['j']) for aln in alignments], [])
        x = np.zeros([len(i), len(j), 5])
        
        # add alignments
        J = 0
        for aln in alignments:
            I = np.array([ind in aln['i'] for ind in i])
            x[I,J:(J+len(aln['j'])),:] = aln['x']
            J = J + len(aln['j'])
        
        # get results
        y = {'x':x, 'i':i, 'j':j}
        with open(os.path.join(args.output_dir, '%s.multigene.aln.pickle' %(genome)), 'wb') as f:
            pickle.dump(y, f)
