import pandas as pd, re, os, sys
from pathlib import Path
import argparse
import glob
import logging

sys.path.append(str(Path(__file__).resolve()))
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

from ut import seq as seq_util
from tools import get_markers

parser = argparse.ArgumentParser()
parser.add_argument('--multigene_table', help='Path to the multigene table')
parser.add_argument('--uhgg_dir', help='Path to the UHGG directory')
parser.add_argument('--uhgg_metadata', help='Name of UHGG metadata', default='genomes-all_metadata.tsv')
parser.add_argument('--genes', metavar='G', type=str, nargs='+',
                        help='List of marker genes')
parser.add_argument('--output_file', default='pangenome_markers.fasta', help='Output file name')
                    
args = parser.parse_args()

n_genes = len(args.genes)

df = pd.read_csv(args.multigene_table, sep='\t')

gi = list(set(df.iloc[:,1:n_genes+1].values.flatten()))
gi = list(set([re.sub('_\d+$','',i) for i in gi]))

x = pd.read_csv(os.path.join(args.uhgg_dir, args.uhgg_metadata), sep='\t',low_memory=False).set_index('Genome')

mg = list(set(x.loc[gi,].MGnify_accession.values))+list(set(df.mgyg))
#mg = list(set(x.MGnify_accession.values))
mg = list(set(mg))


fl = [os.path.join(args.uhgg_dir, 'all_genomes', mi[0:13], mi) for mi in mg]

gff_files = []
for f in fl:
    gff_files.extend(glob.glob(f + "/*/*gff.gz"))

gff_files = list(set(gff_files))
logging.info(f'Found {len(gff_files)} gff files')

# run get_marker_genes.py
get_markers.from_gff(gff_files, args.output_file, args.genes)