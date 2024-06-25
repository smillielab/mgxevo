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
parser$add_argument('--panel')
parser$add_argument('--abd')
parser$add_argument('--strain_freq')
parser$add_argument('--phylo_enrichment')
parser$add_argument('--uhgg_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

uhgg.dir = args$uhgg_dir
out.dir = args$outdir


gi = 'GUT_GENOME147673'
#----------------
# LOAD ABUNDANCES
#----------------

# filtering criteria
cut1 = 0.5
cut2 = 1000

# load bacterial metadata and use on non redundant 
bact = ffread(args$panel, as.dt=T)
setkey(bact, genome)
                                     
# metadata                             
meta = read.metadata.nice(args$metadata)

# gene counts
# -----------
# ab.g = gene counts matrix [5653 samples x 1987 genes]
# note: contains redundant genes
ab.g = ffread(args$abd, row.names=1, header=T)
i = rowSums(ab.g) >= 1e3
ab.g = ab.g[i,]

# refine panel to only include taxa that survive filtering
bact = bact[apply(bact[,c('dnaG','rpoB','gyrB')], 1, function(a) all(a %in% colnames(ab.g))),]
setkey(bact, genome)

# species counts
# --------------
# ab.b = species counts matrix [5653 samples x 243 species]
# note: non-redundant
ab.b = do.call(cbind, sapply(bact$genome, function(gi) rowMeans(ab.g[,as.character(bact[gi,c('dnaG','rpoB','gyrB')])]), simplify=F))
i = rowSums(ab.b) >= 1e4
ab.b = ab.b[i,]


# tpm normalization
# -----------------
tpm.b = 1e4*ab.b/rowSums(ab.b)

# strain counts
# -------------
# ab.s = strain counts matrix [5653 samples x 377 strains]
# note: non-redundant

# load strain frequencies
freq = ffread(args$strain_freq,sep=',',row.names=1, header=T)
freq = freq[, gsub('\\.(.*)$','',colnames(freq)) %in% colnames(ab.b)]

# align matrices
i = intersect(rownames(freq), rownames(ab.b))
freq = freq[i,]
ab.b = ab.b[i,]
ab.s = ab.b[i,]
tpm.b = tpm.b[i,]
tpm.s = tpm.b[i,]

# filter samples
temp = freq * tpm.s[,gsub('\\.(.*)$','',colnames(freq))]
i = rowSums(temp) >= cut2
freq = freq[i,]
ab.b = ab.b[i,]
ab.s = ab.s[i,]
tpm.b = tpm.b[i,]
tpm.s = tpm.s[i,]

f = freq
f[f <= cut1] = 0
freq = f

# strain abundances
ab.s = freq * ab.s[,gsub('\\.(.*)$','',colnames(freq))]
tpm.s = freq * tpm.s[,gsub('\\.(.*)$','',colnames(freq))]
tpm.b = tpm.b[rownames(tpm.s),]

# ----------------
# load strain info
# ----------------

# get tree stats
ph = do.call(rbind, sapply(bact$genome, function(gi) fread(sprintf('%s/%s.phylo.enrich.full.txt', args$phylo_enrichment,gi)) ,simplify=F))
ph = ph[test == 'nice_dxIBD' | is.na(test)]
ph[, padj:=p.adjust(pval, 'fdr'), .(name)]
ph[, id := do.call(paste, c(.SD[,1:2], sep = "."))]
ph = ph[colnames(tpm.s), on=.(id)]
ph[,name:=gsub('\\.(.*)$','', id)]
ph[,node:=gsub('(.*)\\.','', id)]

# get health and disease strains
id.h = ph[padj < .05][coef < 0]$id
id.d = ph[padj < .05][coef > 0]$id
id.b = sort(unique(c(id.h, id.d)))

# --------------------
# aggregate by species
# --------------------

aggbyspecies = function(x, species.use){
    sapply(species.use, function(a){
        j = grep(a, colnames(x))
	if(length(j) > 0){
	    rowSums(x[,j,drop=F])
	} else {
	    setNames(rep(0, nrow(x)), rownames(x))
	}
    })
}

species.use = unique(sort(gsub('\\..*', '', id.b)))

# HC strains
xh = as.matrix(tpm.s[, id.h])
xh = aggbyspecies(xh, species.use)
fh = as.matrix(freq[,id.h])
fh = aggbyspecies(fh, species.use)

# IBD strains
xd = as.matrix(tpm.s[, id.d])
xd = aggbyspecies(xd, species.use)
fd = as.matrix(freq[,id.d])
fd = aggbyspecies(fd, species.use)

# Other strains
xo = as.matrix(tpm.s[, setdiff(colnames(tpm.s), id.b)])
xo = aggbyspecies(xo, species.use)
fo = as.matrix(freq[, setdiff(colnames(tpm.s), id.b)])
fo = aggbyspecies(fo, species.use)

# All strains
xa = as.matrix(tpm.s)
xa = aggbyspecies(xa, species.use)
fa = as.matrix(freq)
fa = aggbyspecies(fa, species.use)

# All species
xb = as.matrix(tpm.b)
                           
# ------------
# cleanup data
# ------------

# log transform
# -------------
i = rownames(xh)
j = colnames(xh)
xh = log2(xh[i,j] + 1)
xd = log2(xd[i,j] + 1)
xo = log2(xo[i,j] + 1)
xa = log2(xa[i,j] + 1)
xb = log2(xb[i,] + 1)
xa[xh + xd + xo == 0] = 0

# fix names
# ---------
colnames(xh) = gsub('^', 'H.', colnames(xh))
colnames(xd) = gsub('^', 'D.', colnames(xd))
colnames(xo) = gsub('^', 'O.', colnames(xo))
colnames(xa) = gsub('^', 'A.', colnames(xa))
colnames(xb) = gsub('^', 'S.', colnames(xb))
                           
# data frame
# ----------
if(paste0('H.', gi) %in% colnames(xh)){
    d = data.frame(
        h = xh[,paste0('H.', gi)],
        d = xd[,paste0('D.', gi)],
        s = xb[,paste0('S.', gi)]
    )
} else {
        d = data.frame(
        h = NA*xb[,paste0('S.', gi)],
        d = NA*xb[,paste0('S.', gi)],
        s = xb[,paste0('S.', gi)]
    )
}
writeTable(d,sprintf('%s/el.ab.tsv', out.dir))
                           
# core genes
# ----------
meta = ffread(paste0(args$uhgg_dir,'/genomes-all_metadata.nice.tsv'))
rownames(meta) = meta$Genome
mi  = meta[gi]$mgyg
ifn = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/genes_presence-absence_locus.mapped.csv',args$uhgg_dir, substr(mi,1,13), mi)
pan = fread(ifn)
pan0 = pan
pan = pan[,-c(1:14)]

pan = sapply(colnames(pan), function(nm) drop.na(unique(unlist(str_split(pan[[nm]], '\t')))), simplify=F)
genes.all = unique(unlist(pan,use.names=F))
U = do.call(cbind, sapply(names(pan), function(nm) genes.all %in% pan[[nm]], simplify=F))
rownames(U) = genes.all
U = Matrix(U, sparse=T)*1
core = data.table(gene=names(rowMeans(U)),pcid=rowMeans(U))
core = core[pcid>.8]
core.genes = unique(core$gene)
writeLines(core.genes,sprintf('%s/el.core.txt', out.dir))