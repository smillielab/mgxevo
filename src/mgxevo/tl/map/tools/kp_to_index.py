import sys

ifn = sys.argv[1]
ofn = sys.argv[2]

x = {}
g = []


flag = 0
for i, line in enumerate(open(ifn)):
    if line.startswith('Sample'):
        flag = 1
        continue
    if flag == 1:
        line = line.rstrip().split('\t')
        gene = line[3]
        if gene not in x:
            x[gene] = [i, i]
            g.append(gene)
        x[gene][1] = i

out = open(ofn, 'w')
for gene in g:
    out.write('%s\t%s\t%s\n' %(gene, x[gene][0], x[gene][1]))
out.close()
