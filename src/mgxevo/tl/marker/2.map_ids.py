gmap = {}
import argparse

def process_genes(genes, uclust_path, output_file):
    gmap = {}

    for gene in genes:
        for line in open(f'{uclust_path}/{gene}.uc'):
            line = line.rstrip().split('\t')
            hit = line[8]
            seed  = line[9]
            if seed == '*':
                seed = hit
            gid = '_'.join(hit.split('_')[:2])
            if gid not in gmap:
                gmap[gid] = {}
            if gene not in gmap[gid]:
                gmap[gid][gene] = []
            gmap[gid][gene].append(seed)

    with open(output_file, 'w') as f:
        f.write('\t'.join(['genome'] + genes) + '\n')
        for genome in gmap:
            out = [genome]
            for gene in genes:
                gids = ','.join(list(set(gmap[genome].get(gene, []))))
                out.append(gids)
            f.write('\t'.join(out) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process marker genes.')
    parser.add_argument('--uclust_path', type=str, default='uclust',)
    parser.add_argument('--genes', metavar='G', type=str, nargs='+',
                        help='a list of marker genes to be processed')
    parser.add_argument('--output_file', type=str, default='gene_ids.txt',
                        help='the output file where the results will be written')
    args = parser.parse_args()
    process_genes(args.genes, args.uclust_path, args.output_file)