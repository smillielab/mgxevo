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
        
    library(caret)
    library(randomForest)
    library(gbm)

    })
})

parser = ArgumentParser()
parser$add_argument("--metadata")
parser$add_argument('--panel')
parser$add_argument('--nr')
parser$add_argument('--abd')
parser$add_argument('--freqs')
parser$add_argument('--phylo_enrichment')
parser$add_argument('--outdir')
parser$add_argument('--cut1')
parser$add_argument('--cut2')
parser$add_argument('--ibd')
parser$add_argument('--fopt')
parser$add_argument('--coh')
parser$add_argument('--out')
args = parser$parse_args()

run_m1 = function(x, i=NULL){

    if(is.null(i)){
        i = rownames(x)
    }
    
    # model 1: hc vs ibd
    x1 = x[i,]
    y1 = factor(ifelse(meta[i,]$nice_dx2 == 'HC', 'HC', 'IBD'), levels=c('HC', 'IBD'))
    r1 = randomForest(x1, y1, sampsize=c(500,500), importance=T)
    R1 = de.fdr(t(x1), y1, sens_cut=0)
    
    # model 2: cd vs uc
    x2 = x[i,]
    y2 = factor(meta[i,]$nice_dx2, levels=c('CD', 'UC'))
    ii = !is.na(y2)
    r2 = randomForest(x2[ii,], y2[ii], sampsize=c(200,200), importance=T)
    R2 = de.fdr(t(x2[ii,]), y2[ii], sens_cut=0)
    
    list(r1=r1, r2=r2, R1=R1, R2=R2)
}

run_m2 = function(x, i=NULL){

    if(is.null(i)){
        i = rownames(x)
    }
    
    # model 1: hc vs ibd
    y1 = factor(ifelse(meta[i,]$nice_dx2 == 'HC', 'HC', 'IBD'), levels=c('HC', 'IBD'))
    
    # model 3: fecal calprotectin (rf)
    x3 = x[i,]
    y3 = log2(fcal[i,]$fcal+1)
    if(ibd){
        ii = (y1 == 'IBD') & (!is.na(y3))
    } else {
        ii = (!is.na(y3))
    }
    print(dim(x3[ii,]))
    r3 = randomForest(x3[ii,], scale(y3[ii]), importance=T)
    R3 = cor(x3[ii,], scale(y3[ii]), method='spearman')

    # model 4: fecal calprotectin (gbm)
    d = data.frame(x3[ii,], y=scale(y3[ii]))
    r4 = gbm(y ~ ., distribution='gaussian', data=d, n.tree=1e3, cv.folds=5, verbose=T)
    #r4 = list()

    list(r3=r3, R3=R3, r4=r4)
}

cut1 = as.numeric(args$cut1)
cut2 = as.numeric(args$cut2)
ibd  = as.logical(args$ibd )
fopt = as.integer(args$fopt)
coh  = args$coh 
out  = args$out 

meta = read.metadata.nice(args$metadata)
rownames(meta) = meta$sample
i = rownames(meta)

fcal = meta


# filter by cohort
if(coh == 'prism'){
    fcal[grep('PRISM', fcal$cohort, invert=T),]$fcal = NA
} else {
    quit('invalid cohort')
}

# load bacterial metadata and use on non redundant 
bact = ffread(args$panel, as.dt=T)
bact = bact[genome %in% readLines(args$nr)]

# gene counts
# -----------
# ab.g = gene counts matrix [5653 samples x 1987 genes]
# note: contains redundant genes
ab.g = ffread(args$abd, sep=',', row.names=1, header=T)
i = rowSums(ab.g) >= 1e3
ab.g = ab.g[i,]
print('ab.g')
print(dim(ab.g))

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
print('ab.b')
print(dim(ab.b))

# tpm normalization
# -----------------
tpm.b = 1e4*ab.b/rowSums(ab.b)

# strain counts
# -------------
# ab.s = strain counts matrix [5653 samples x 377 strains]
# note: non-redundant

# load strain frequencies
freq = ffread(args$freq,sep=',',row.names=1, header=T)
freq = freq[, gsub('\\.(.*)$','',colnames(freq)) %in% colnames(ab.b)]

# align matrices
i = intersect(rownames(freq), rownames(ab.b))
freq = freq[i,]
ab.b = ab.b[i,]
ab.s = ab.b[i,]
tpm.b = tpm.b[i,]
tpm.s = tpm.b[i,]
print('align')
print(dim(freq))

# filter samples
temp = freq * tpm.s[,gsub('\\.(.*)$','',colnames(freq))]
i = rowSums(temp) >= cut2
freq = freq[i,]
ab.b = ab.b[i,]
ab.s = ab.s[i,]
tpm.b = tpm.b[i,]
tpm.s = tpm.s[i,]
print('filter')
print(dim(freq))
                             
# option 1 = filter strain frequencies (no aggregation)
# -----------------------------------------------------

if(fopt == 1){
    f = freq
    f[f <= cut1] = 0
    #g = gsub('\\..*', '', colnames(freq))
    #f = sapply(unique(g), function(a) t(apply(freq[, g == a], 1, function(b){b[b != max(b)] = 0; b})))
    #f = do.call(cbind, f)[,colnames(freq)]
    freq = f
}

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

# option 2 = filter strain frequencies (aggregating by species)
# -------------------------------------------------------------
if(fopt == 2){
    xh[fh <= cut1] = 0
    xd[fd <= cut1] = 0
    xo[fo <= cut1] = 0
}

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

# -------------------
# disease predictions
# -------------------

# create input data
# -----------------

x = list()
x$h = xh
x$d = xd
x$o = xo
x$a = xa
x$b = xb
x$hd = cbind(xh, xd)
x$ha = cbind(xh, xa)
x$da = cbind(xd, xa)
x$oa = cbind(xo, xa)
x$hdo = cbind(xh, xd, xo)
x$hda = cbind(xh, xd, xa)
x$all = cbind(xh, xd, xo, xa)

# run predictions
# ---------------

g = meta[rownames(x$h),]$subject
i = tapply(rownames(x$h), g, sample, 1)
m1 = sapply(x, function(a) run_m1(a, i=i), simplify=F)

g = meta[rownames(x$h),]$subject
j = is.na(fcal[rownames(x$h),]$fcal)
g[j] = NA
i = tapply(rownames(x$h), g, sample, 1)
print('fcal')
print(length(i))
m2 = sapply(x, function(a) run_m2(a, i=i), simplify=F)

saveRDS(list(m1=m1, m2=m2, i=i), file=paste0(args$outdir,'/',args$out))