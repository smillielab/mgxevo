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
parser$add_argument("--abd")
parser$add_argument("--uhgp_tools")
parser$add_argument('--dmd')
parser$add_argument('--phylo_enrichment')
parser$add_argument('--uhgg_dir')
parser$add_argument('--tree_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

# get tree stats
gi = el
ph = fread(sprintf('%s/%s.phylo.enrich.full.txt', args$phylo_enrichment,gi))
ph = ph[test == 'nice_dxIBD' | is.na(test)]
ph[, padj:=p.adjust(pval, 'fdr'), .(name)]
ph[, id := do.call(paste, c(.SD[,1:2], sep = "."))]
ph[,name:=gsub('\\.(.*)$','', id)]
ph[,node:=gsub('(.*)\\.','', id)]
setkey(ph, id)

# get health and disease strains
id.h = ph[padj < .05][coef < 0]$id
id.d = ph[padj < .05][coef > 0]$id

panel = data.table(expand.grid(h=id.h, d=id.d))[,h.pval := ph[h]$padj][,d.pval:=ph[d]$padj]
panel = panel[,.SD[which.min(d.pval*h.pval)]]
panel[,h:=as.numeric(gsub('(.*)\\.','', h))][,d:=as.numeric(gsub('(.*)\\.','', d))]

# get samples
tr = read.tree.nice(sprintf('%s/%s.tree',args$tree_dir, gi))
h = strip.dx_labs(extract.clade(tr, node = panel$h)$tip.label)
d = strip.dx_labs(extract.clade(tr, node = panel$d)$tip.label)
i = unique(c(h,d))

# load matrix
X  = readRDS(arg$dmd)
X = Matrix::Diagonal(x = 1 / Matrix::rowSums(X)) %*% X
X@x = log2((1e4 * (X@x)) + 1)


# approach 2: add E. lenta abundance as a covariate
# -------------------------------------------------
Q = X

ab = ffread(args$abd)
i = intersect(i, rownames(ab))

Q = Q[i,]
dg = data.frame(sample=i,dx=factor(ifelse(i %in% h, 'HC', 'IBD'),levels=c('HC', 'IBD')))

Q = as.matrix(Q)
Q = Q[dg$sample,]
Q = scale(Q)
Q = Q[,which.names(!apply(is.na(Q),2, all))]

dx = dg$dx 
res = lm(Q[dg$sample,] ~ ab[dg$sample,'s'] + dx)
res = summary(res)
res = sapply(names(res), function(nm) cbind(gene = gsub('(.*)\\s','',nm), 
                                            res[[nm]]$coef['dxIBD',c('Estimate', 'Pr(>|t|)'),drop=F]), simplify=F)
res = data.table(do.call(rbind, res))
colnames(res) = c('gene', 'coeff', 'pval')
res[,coeff:=as.numeric(coeff)]
res[,pval:=as.numeric(pval)]
res[,padj:=p.adjust(pval,'fdr')]

x2 = res
             
#--------------------
# Map onto DE results
#--------------------
                      
# map UHGP90, KEGG-KOs and annotation groups
A = as.data.table(stack(kegg$g))[,.(ko=values, anno=ind)]
B = fread(sprintf('%s/uhgp-90.ko', args$uhgg_dir))
B = merge(A, B, by='ko', allow.cartesian = T)
B = B[,.(gene=id, ko, anno)]

# table instances of each group
B = B[gene %in% res$gene]
B = B[,.(.N, ko),.(gene, anno)]
B = merge(B, A[,.(N2=.N),.(anno)], by='anno', all.x=T)
B = B[,.(gene, ko, anno, N, N2)]
B = merge(B[,.(ko=paste(unique(ko), collapse=', ')),.(gene)], B[,.SD[order(-N, N2)][1], .(gene)][,.(gene, anno)], by='gene')

# merge
x2 = merge(x2, B, by='gene', all.x=T)
cmd = sprintf('python %s/tools/clean_eggNOG.py %s/uhgp-90_eggNOG.tsv', args$uhgp_tools, args$uhgg_dir)
gene.symbol = fread(cmd=cmd,select = c(1,2))
gene.symbol[name=='',name:=NA]
x2 = merge(x2, gene.symbol, by='gene', all.x=T)

fwrite(x2, sprintf('%s/elenta.de.csv', args$outdir))