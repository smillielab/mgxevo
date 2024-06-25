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
parser$add_argument("--output_dir", help="Output directory path")
args <- parser$parse_args()

abundance_file <- args$abundance_file
metadata <- args$metadata
map_ids <- args$map_ids
uhgg_catalogue <- args$uhgg_catalogue
uhgg_metadata <- args$uhgg_metadata
genes <- args$genes

# load data
x = ffread(abundance_file, row.names=1, header=T, sep=',')
dx = read.metadata.nice(metadata)
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
set.seed(3)
noise = rnorm(n=nrow(x)*ncol(x), mean=1e-5, sd=1e-5)
x = x + abs(noise)

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
    u = sapply(m[i,], function(a) sample(a, size=1), simplify=F)
    U = as.character(unlist(u))
    V = colMeans(x[,U])
    W = colSums(x[,U] >= 10)
    c(U, V, W)

}, simplify=F)

# convert to dataframe
QQ = data.frame(do.call(rbind, Q))
rownames(QQ) = names(Q)
    
u_genes = paste0('u_', genes)
n_genes = paste0('n_', genes)
colnames(QQ) = c(genes, u_genes, n_genes)
    
# fix numeric
QQ[,(length(genes)+1):ncol(QQ)] = sapply(QQ[,(length(genes)+1):ncol(QQ)], as.numeric)
QQ = as.data.table(QQ %>% rownames_to_column('genome'))

# address duplicate markers
# -------------------------
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

# remove duplicates, selecting the genomes with the greatest number of annotations 
x = x[order(-x$num.anno),]
x = x[!duplicated(x[,genes]),]

# append useful info
# ------------------
    
# add "nice" names and extra species metadata 
v = x
v$ab = rowMeans(x[,u_genes])
m = ffread(uhgg_metadata)
anno = gsub('^(d|p|c|o|f|g|s)__','',str_split_fixed(m[match(v$genome,m$Genome),]$Lineage, ';', n=7))
colnames(anno) = c('domain','phylum','class','order', 'family','genus', 'species')

v = cbind(v,anno)
                                                
j = c('genome', genes, u_genes, n_genes,'ab','domain','phylum','class','order','family','genus','species','mgyg','num.anno')
v = v[,j]
fwrite(v, paste0(args$output_dir, '/ids.multigene.big.txt'), row.names=F, sep='\t', quote=F)
