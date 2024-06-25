require(data.table)
require(optparse)

option_list = list(make_option('--name', help='name', default='test'),
	           make_option('--dge', help='dge filename'),
	           make_option('--sep', help='field separator', default='\\.'),
	           make_option('--field', help='collection field', type='integer', default=1)
)

args = parse_args(OptionParser(option_list=option_list))

# Get filenames
dge = data.frame(fread(paste('zcat', args$dge)), row.names=1)

# Get collections
collections = sapply(strsplit(colnames(dge), args$sep), '[', args$field)
print(paste('Found', length(unique(collections)), 'collections'))
u = unique(collections)

# Calculate quality statistics
quals1 = c() # per-collection
quals2 = c() # per-experiment
for(cutoff in c(500, 750, 1000)){

    qual1 = matrix(0, nrow=length(u), ncol=3)
    qual2 = rep(0, 4)
    rownames(qual1) = u

    num_genes = colSums(dge > 0)
    num_reads = colSums(dge)
    cells = tapply(num_genes >= cutoff, collections, sum)
    genes = tapply(num_genes, collections, median)
    reads = tapply(num_reads, collections, median)
    qual1[u, 1:3] = cbind(cells[u], genes[u], reads[u])
    qual2[1:4] = c(sum(num_genes >= cutoff), sum(num_genes >= cutoff)/length(u), median(num_genes), median(num_reads))

    quals1 = cbind(quals1, qual1)
    quals2 = c(quals2, qual2)
}

q1 = expand.grid(c('Total cells', 'Median genes per cell', 'Median reads per cell'), c('(cutoff=500)', '(cutoff=750)', '(cutoff=1000)'))
q1 = c('Project', paste(as.character(q1[,1]), as.character(q1[,2])))

q2 = expand.grid(c('Total cells', 'Cells per collection', 'Median genes per cell', 'Median reads per cell'), c('(cutoff=500)', '(cutoff=750)', '(cutoff=1000)'))
q2 = c('Project', paste(as.character(q2[,1]), as.character(q2[,2])))

quals1 = cbind(rep(as.character(args$name), nrow(quals1)), quals1)
colnames(quals1) = q1

quals2 = c(args$name, quals2)
quals2 = t(data.frame(cbind(q2, quals2), row.names=1))

write.table(quals1, file=paste(args$name, '.collections.dge_qual.txt', sep=''), sep='\t', quote=F)
write.table(quals2, file=paste(args$name, '.experiments.dge_qual.txt', sep=''), sep='\t', quote=F, row.names=F)
