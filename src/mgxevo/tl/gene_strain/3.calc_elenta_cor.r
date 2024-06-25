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
parser$add_argument("--metadata")
parser$add_argument('--el_abd')
parser$add_argument('--dmd_dir')
parser$add_argument("--gene100", type="character", help="Genes over 100 counts file")
parser$add_argument('--outdir')
parser$add_argument("--chunk", type="integer", help="Index")

args = parser$parse_args()


genes = sort(unique(readLines(args$gene100)))

m = read.metadata.nice(args$metadata)
rownames(m) = m$sample

x = ffread(args$el_abd)
x = Matrix(as.matrix(x),sparse=T)
x = pad_mtx(x, rows = m$id, cols = colnames(x))

chunk2 = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
col.use = chunk2(genes, 50)[[args$chunk]]

rows = rownames(x)
cols = col.use

M = Matrix(0, nrow=length(rows), ncol=length(cols), sparse=TRUE)
rownames(M) = rows
colnames(M) = cols

# iterate through DIAMOND files
for(n in 1:50){
    M =  M + pad_mtx(readRDS(sprintf('%s/diamond_counts.i%s.n50.rds',args$dmd_dir, n)), rows = rows, cols = cols)
}

ROW_SUMS = readRDS(sprintf('%s/rowSums.dmd.rds', args$outdir))

i = which.names(ROW_SUMS >= 1e4)
M = M[i,]

i = intersect(rownames(M), rownames(x))
M = M[i,]
x = x[i,]

# tpm for genes
g = M/ROW_SUMS[rownames(M)]
g@x = log2(1e4*g@x + 1)

# ------------          
# correlations
# ------------
            
# mask 0-0 pairs for correlations with NA
mask = as(!sweep_sparse(g, 1, x[,'s'], '+'),'lgCMatrix')
idx = which(mask)

g.mask = g
g.mask[idx] = NA
rho = cor(x[,'s'], g.mask, use='pairwise.complete.obs', method='spearman')
ndf = colSums(!is.na(g.mask))
p.rho = 1-pf(rho^2*(ndf-2)/(1- rho^2), 1, ndf-2)
i = colnames(rho)
cor.res = data.table(gene = i, rho = rho[1,i], pval = p.rho[1,i], n.obs = ndf[i])
out = sprintf('%s/elenta.cor.i%s.tsv',args$outdir,args$chunk)
fwrite(cor.res,out, row.names=FALSE)

