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
parser$add_argument('--gene_strain_ids')
parser$add_argument('--uhgg_dir')
parser$add_argument('--mtx_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

gi = args$genome
uhgg.dir = args$uhgg_dir
ids.use.fn = args$gene_strain_ids
mtx.dir = args$mtx_dir
out.dir = args$outdir

#--------------
# load metadata
#--------------
meta = ffread(sprintf('%s/genomes-all_metadata.nice.tsv',uhgg.dir))
rownames(meta) = meta$Genome

ids = fread(ids.use.fn)

# read data
ifn = sprintf('%s/%s.mtx.rds',mtx.dir, gi)
dat = readRDS(ifn)

# mapping
mp = data.table(dat[['d']])
key = list(K='kegg', G='gene', U = 'uhgp')

compare.against= 'D'
if(ids[genome==gi,test]=='other(D)/H') compare.against = 'H'

res = sapply(c('K','G', 'U'), function(i){
    M = as(dat[[i]], "lMatrix")
    # fisher test for IBD vs HC (odds_ratio > 1 ~ IBD)
    v = factor(mp[colnames(M), on=.(gid)]$type==compare.against, levels=c(FALSE,TRUE))
    
    q = sapply(1:nrow(M), function(i){
        u = factor(M[i,], levels=c(FALSE,TRUE))
        ft = fisher.test(table(u,v))
        c(odds_ratio=unname(ft$estimate), pval = ft$p.value)} ,simplify=F)
    q = cbind(gene=rownames(M), data.frame(do.call(rbind,q)), matrix=key[[i]])
    
    # prob of IBD and HC
    p = sapply(unique(mp$type), function(vi) rowMeans(M[,mp[mp$type==vi,]$gid,drop=F]))
    p = data.frame(p)
    colnames(p) = paste0('prob.',tolower(substr(colnames(p),1,1)))
    p = cbind(p,gene=rownames(p))
    
    # merge
    q = data.table(merge(q, p, by='gene',all=T))
    q[,padj:=p.adjust(pval,'fdr')]
    q[,genome:=gi]
    q},simplify=F)
res = do.call(rbind,res)
    
# write output
ofn = sprintf('%s/%s.txt', out.dir, gi)
fwrite(res, ofn, sep='\t', quote=F, row.names = F)