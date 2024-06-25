import argparse, os, re, sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

import ut

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gene', help='gene name')
parser.add_argument('--multi', help='print multiple genes?', default=False, action='store_true')
parser.add_argument('--mgyg', help='mgyg list (core only)', default='uhgg/genomes-nr.mgyg_ids.txt')
parser.add_argument('--regex', help='mgyg regex (core only)', default='')
parser.add_argument('--uhgg', help='uhgg catalogue', default='uhgg')
parser.add_argument('--out', help='output fasta', default='uhgg/marker_genes/genes.fna')
args = parser.parse_args()

args.out = args.out.replace('genes.fna', '%s.fna' %(args.gene))

# get fasta sequences from core genome

# get genome ids
# --------------
# option 1: mygyg list
gids = [line.rstrip() for line in open(args.mgyg)]
# option 2: mgyg regex
if args.regex:
    gids = []
    for line in open(f'{args.uhgg}/genomes-nr_metadata.tsv'):
        [mgyg, lineage] = line.rstrip().split('\t')[17:19]
        m = re.search(args.regex, lineage)
        if m:
            gids.append(mgyg)
    gids = sorted(list(set(gids)))
    logging.info('Found %d mgyg ids using regex %s:\n%s' %(len(gids), args.regex, ' '.join(gids)))
    
    
# get fasta sequences
os.makedirs(os.path.dirname(args.out), exist_ok=True)

with open(args.out, 'w') as outfile:
    for gid in gids:

        # get filenames
        pre = gid[:13]
        gff = f'{args.uhgg}/uhgg_catalogue/%s/%s/genome/%s.gff' %(pre, gid, gid)
        fna = f'{args.uhgg}/uhgg_catalogue/%s/%s/genome/%s.fna' %(pre, gid, gid)
        
        # check if files exist
        if not os.path.exists(gff) or not os.path.exists(fna):
            sys.stderr.write('error: could not find %s\n' %(fna))
            continue
        
        # store gene coordinates
        # - coords[contig] = [beg, end, strand, gene]
        coords = {}
        
        # find gene coordinates in gff file
        for line in open(gff):
            m = re.search('gene=%s[;_]' %(args.gene), line)
            if m:

                # extract gene coordinates
                line = line.rstrip().split('\t')
                contig = line[0].strip()
                beg = int(line[3])
                end = int(line[4])
                strand = line[6]
                gene = line[8].split(';')[0][3:].strip()
                
                # store gene coordinates
                if contig not in coords:
                    coords[contig] = []
                coords[contig].append([beg, end, strand, gene])
                
                # store one coord unless --multi
                if not args.multi:
                    break
        
        # skip empty results
        if len(coords) == 0:
            continue
            
        # write fasta sequences
        for record in ut.iter_fst(fna):
            contig = record[0][1:]
            if contig in coords:
                
                # write each gene sequence
                for coord in coords[contig]:
                    [beg, end, strand, gene] = coord
                    seq = record[1][(beg-1):end]
                    if strand == '-':
                        seq = ut.reverse_complement(seq)
                    outfile.write('>%s\n%s\n' %(gene, seq))
                
                # write one sequence unless --multi
                if not args.multi:
                    break