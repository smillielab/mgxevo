import argparse
import glob
import gzip
import re
import logging

def open_file(file_name):
    if isinstance(file_name, str):
        if '.gz' in file_name:
            file_name = gzip.open(file_name, 'rt')
        else:
            file_name = open(file_name)
    return file_name

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', help='BLAST file', required=True)
    parser.add_argument('--pct', help='Percent identity', default=90, type=float)
    parser.add_argument('--evalue', help='Evalue cutoff', default=1e-5, type=float)
    parser.add_argument('--output', help='Output prefix (counts)')
    return parser.parse_args()

def parse_blast_report(args, count):
    logging.info('Parsing BLAST report')
    for line in open_file(args.blast):
        line = line.rstrip().split('\t')
        target = line[1]
        pct = float(line[2])
        evalue = float(line[10])
        if pct <= args.pct or evalue >= args.evalue:
            continue
        count[target] = count.get(target, 0) + 1
    return count

def write_otu_tables(args, count):
    logging.info('Writing OTU tables')
    with open('%s.counts.txt' %(args.output), 'w') as out:
        oi = ['Gene'] + sorted(count)
        out.write('\t'.join(oi) + '\n')
        oi = [args.output] + [str(count[k]) for k in sorted(count)]
        out.write('\t'.join(oi) + '\n')
        
def main():
    args = parse_arguments()
    count = {}
    count = parse_blast_report(args, count)
    write_otu_tables(args, count)

if __name__ == "__main__":
    main()
