script_path <- this.path::this.path()
script_dir <- dirname(script_path)

# SELECT IDS FOR MULTIGENE PANEL
util <- file.path(script_dir, "/../../ut/util.r")
treetools <- file.path(script_dir, "/../tree/treetools.r")
source(util)
source(treetools)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--abundance_file", help="Abundance file path")
parser$add_argument("--metadata", help="Metadata file path")
parser$add_argument("--map_ids", help="Map IDs file path")
parser$add_argument("--uhgg_catalogue", help="UHGG catalogue path")
parser$add_argument("--uhgg_metadata", help="UHGG metadata file path")
parser$add_argument("--genes", type="character", nargs='+',
                        help='a list of marker genes to be processed')
parser$add_argument("--rho_cutoff", type="numeric", default=0.5, help="Rho cutoff value")
parser$add_argument("--nsample_cutoff", type="numeric", default=50, help="Minimum number of samples with 10 or more read counts")
parser$add_argument("--output_dir", help="Output directory path")
args <- parser$parse_args()

abundance_file <- args$abundance_file
metadata <- args$metadata
map_ids <- args$map_ids
uhgg_catalogue <- args$uhgg_catalogue
uhgg_metadata <- args$uhgg_metadata
genes <- args$genes
nsample_cutoff <- args$nsample_cutoff

# load data
x = ffread(abundance_file, row.names=1, header=T, sep=',')
dx = read.metadata.nice(args$metadata)
rownames(dx) = dx$sample

if(nrow(dx) != nrow(x)) {
  warning("Number of rows in dx and x do not match.")
}
dx = dx[rownames(x), ]

# aggregate by samples by subject
x = nice_agg(as.matrix(x), dx$subject, type = 'mean')

# exclude subjects that have fewer than 10 counts across all genes
i = rowSums(x) > 10
x = x[i,]

# only genes where at least 10 samples have 10 or more read counts
exlcude = colSums(x >= 10) <= 10
x = x[,!exlcude]
exlcude = which.names(exlcude)

# add noise to fix high correlations with low abundance samples
noise = rnorm(n=nrow(x)*ncol(x), mean=1e-5, sd=1e-5)
x = x + abs(noise)

# TP10K
# x = log2(1e4*x/rowSums(x) + 1)
x = 1e4*x/rowSums(x)

# from map_ids.py
m = ffread(map_ids, row.names=1, header=T)

# fix mapping file to include genes that intersect with `x`
m = data.frame(apply(m, c(1,2), function(a) paste0(setdiff(unlist(strsplit(a,',')), exlcude),collapse=',')))

# remove empty
i = apply(m != '', 1, all)
m = m[i,]

# split genes
i = rownames(m)
m = sapply(m, strsplit, ',')
rownames(m) = i

# select genes
Q = sapply(rownames(m), function(i){
    # enumerate possibilities
    u = expand.grid(m[i,])

    # calculate all correlations
    v = intersect(sort(unique(unlist(m[i,]))), colnames(x))
    q = na.replace(cor(x[,v], method='spearman'), 0)
    
    # find maximum correlation
    z = apply(u, 1, function(a) mean(q[a,a]))
    U = as.character(unlist(u[which.max(z),]))
    V = colMeans(x[,U])
    W = colSums(x[,U] >= 10)
    c(U, V, W, max(z))

}, simplify=F)

# convert to dataframe
QQ = data.frame(do.call(rbind, Q))
rownames(QQ) = names(Q)

u_genes = paste0('u_', genes)
n_genes = paste0('n_', genes)
colnames(QQ) = c(genes, u_genes, n_genes, 'rho')

# fix numeric
QQ[,4:ncol(QQ)] = sapply(QQ[,4:ncol(QQ)], as.numeric)
QQ = as.data.table(QQ %>% rownames_to_column('genome'))
    
x = data.frame(QQ)
m = ffread(uhgg_metadata, sep='\t', row.names=NULL)
# remove genomes that have identical panels
x$mgyg = m[match(x$genome,m$Genome),]$MGnify_accession

num.anno = sapply(x$mgyg, function(mg){
    df = ffread(paste0(uhgg_catalogue, '/', str_sub(mg,1,13),'/',mg,'/genome/',mg,'.anno.nice.txt'))
    df = df[,c('eggnog_desc', 'interpo_desc', 'accesion_desc')]
    df[df=='']=NA
    df = df[complete.cases(df),]
    i = !apply(sapply(colnames(df), function(a) grepl('unknown|hypothetical|uncharacterized', df[,a],ignore.case=TRUE)),1,any)
    df = df[i,]
    nrow(df)
},simplify =T)
x$num.anno = num.anno[x$mgyg] 
m = ffread(uhgg_metadata)
count = nrow(x)

# remove duplicates, selecting the genomes with the greatest number of annotations 
x = x[order(-x$num.anno),]
x = x[!duplicated(x[,genes]),]
count = nrow(x)

# keep only well correlated genomes
rho.cutoff = args$rho_cutoff
x = x[x$rho > rho.cutoff,]
count = nrow(x)

v = x
v$ab = rowMeans(x[,u_genes])

anno = gsub('^(d|p|c|o|f|g|s)__','',str_split_fixed(m[match(v$genome,m$Genome),]$Lineage, ';', n=7))
colnames(anno) = c('domain','phylum','class','order', 'family','genus', 'species')

v = cbind(v,anno)
                                                
j = c('genome', genes, u_genes, n_genes, 'rho','ab','filter.min.ab','domain','phylum','class','order','family','genus','species','mgyg','num.anno')
v = v[,j]

# nice names
nm = v[,c('domain','phylum','class','order', 'family','genus', 'species')]
nm[nm$species=='',]$species = paste0(nm[nm$species=='',]$genus, ' sp.')
if(nrow(nm[nm$species==' sp.',]) > 0){
    nm[nm$species==' sp.',]$species = ''
}
nice.name = sapply(nm, function(a) !grepl("\\d",a) & !(a==''))
ix = sapply(1:nrow(nm), function(a) max(which(nice.name[a,])))

nice.name = sapply(1:nrow(nm),function(a){
    if(ix[a]>5){
        if(ix[a]==6){
            paste0(nm[a,]$genus, ' sp.')
        } else {
            nm[a,]$species
        }
    } else {
        if(nm[a,]$species!=''){
            paste0(nm[a,]$species,' (',nm[a,ix[a]],' ',colnames(nm)[ix[a]],')')
        } else{
            paste0(nm[a,ix[a]],' ',colnames(nm)[ix[a]],' microbe')
        }
    }
})
v$nice_name = nice.name

# remove the 2 at the end of the txt if not using log2(TP10K+1)
fwrite(v, paste0(args$output_dir, '/ids.multigene.txt'), row.names=F, sep='\t', quote=F)
