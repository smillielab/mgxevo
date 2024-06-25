suppressWarnings({
    suppressMessages({
    library(argparse)

    script_path <- this.path::this.path()
    script_dir <- dirname(script_path)
        
    mtx <- file.path(script_dir, "/../../ut/mtx.r")
    treetool <- paste0(script_dir,'/../tree/treetools.r')
    util <- paste0(script_dir,'/../../ut/util.r')
        
    source(mtx)
    library(data.table)
    source(treetool)
    source(util)
    })
})

parser = ArgumentParser()
parser$add_argument('--genome')
parser$add_argument('--metadata')
parser$add_argument('--gene_strain_ids')
parser$add_argument('--uhgg_dir')
parser$add_argument('--phylo_dir')
parser$add_argument('--dmd_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

gi = args$genome
uhgg.dir = args$uhgg_dir
ids.use.fn = args$gene_strain_ids
dmd.dir = args$dmd_dir
out.dir = args$outdir
tree.dir = args$phylo_dir

#--------------
# load metadata
#--------------
ids = ffread(ids.use.fn)
rownames(ids) = ids$genome

meta = ffread(sprintf('%s/genomes-all_metadata.tsv',uhgg.dir))
rownames(meta) = meta$Genome

#----------
# get genes
#----------
# read tree 
ifn = sprintf('%s/%s.ref.tree',gi)
tree = read.tree.nice(ifn)

# get reference genomes
gx = unique(c(gi, fix_ref_names(grep('^GUT', tree$tip.label,value=T))))
mg = unique(meta[gx,]$mgyg)

# get all genes
genes = sapply(mg, function(mi){
    ifn = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/genes_presence-absence_locus.mapped.csv', uhgg.dir, substr(mi,1,13), mi)
    if(file.exists(ifn)){
        genomes = colnames(fread(ifn, header=T, nrow=0))
        j = sort(match(intersect(gx,genomes), genomes))
        mm = bigreadr::big_fread2(ifn, select = j)
        mm  = sapply(mm, function(a) unique(unlist(strsplit(a,'\t'),use.names=F)),simplify=F)
    } else {
        ifn2 = sprintf('%s/uhgg_catalogue/%s/%s/genome/%s.uhgp.txt', uhgg.dir, substr(mi,1,13),mi, mi)
        mm = readLines(ifn2)
        mm = setNames(list(mm),meta[meta$mgyg==mi,]$Genome)
    }
    mm},simplify=F)
genes = unique(unlist(genes,use.names=F))

#-------------
# build matrix
#-------------
m = read.metadata.nice(args$metadata)
rows = m$id
cols = genes

# initialise matrix
M = Matrix(0, nrow=length(rows), ncol=length(cols), sparse=TRUE)
rownames(M) = rows
colnames(M) = cols

# iterate through DIAMOND files
for(n in 1:50){
    M =  M + pad_mtx(readRDS(sprintf('%s/diamond_counts.i%s.n50.rds',dmd.dir,n)), rows = rows, cols = cols)
}
saveRDS(M, file = sprintf('%s/%s.counts.rds',out.dir,gi))