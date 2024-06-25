import argparse, sys, glob, os, re
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import concurrent.futures
import tempfile
import shutil

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from ut import seq as seq_util

def iter_gff(fn, split=False):
    sid = ''
    seq = ''
    i = 0
    for line in seq_util.open_gz(fn):
        i += 1
        if line.startswith('>'):
            i = 1
        line = line.rstrip()
        if line.startswith('>'):
            if (seq != '') & (sid != ''):
                yield [sid, seq]
            sid = line
            if split is True:
                sid = sid.split()[0]
            seq = ''
        else:
            seq += line
    yield [sid, seq]


def process_gff(ifn, marker):
    coords = {}
    for line in seq_util.open_gz(ifn):
        m = any([re.search('gene=%s[;_]' %(gi), line) for gi in marker])
        if m:
            line = line.rstrip().split('\t')
            contig = line[0].strip()
            beg = int(line[3])
            end = int(line[4])
            strand = line[6]
            gene = line[8].split(';')[0][3:].strip()

            if contig not in coords:
                coords[contig] = []
            coords[contig].append([beg, end, strand, gene])
    if len(coords) == 0:
        return

    temp_file = tempfile.NamedTemporaryFile(delete=False)
    for record in iter_gff(ifn):
        contig = record[0][1:]
        if contig in coords:
            for coord in coords[contig]:
                [beg, end, strand, gene] = coord
                seq = record[1][(beg-1):end]
                if strand == '-':
                    seq = seq_util.reverse_complement(seq)
                with open(temp_file.name, 'a') as ofh:
                    ofh.write('>%s\n%s\n' %(gene, seq))
    return temp_file.name

def from_gff(gff_files, output_file, marker=['rpoB','dnaG','gyrB'], threads=None):
    assert isinstance(gff_files, list), '`gff_files` must be a list'
    
    if threads is None:
        threads = os.cpu_count() -2
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_gff, ifn, marker) for ifn in gff_files]
        temp_files = []
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            temp_files.append(future.result())
        
        with open(output_file, 'w') as ofh:
            for temp_file in temp_files:
                if temp_file is None:
                    continue
                
                with open(temp_file, 'r') as ifh:
                    shutil.copyfileobj(ifh, ofh)
                os.remove(temp_file)
            
                
                        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract marker genes from gff files')
    parser.add_argument('-i', '--input', help='input gff files', nargs='+', required=True)
    parser.add_argument('-o', '--output', help='output fasta file', required=True)
    parser.add_argument('-m', '--marker', help='marker genes', nargs='+', default=['rpoB','dnaG','gyrB'])
    args = parser.parse_args()
    from_gff(args.input, args.output, args.marker)