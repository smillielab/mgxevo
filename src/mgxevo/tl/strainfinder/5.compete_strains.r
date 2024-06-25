script_path <- this.path::this.path()
script_dir <- dirname(script_path)

util <- paste0(script_dir,'/../../ut/util.r')
treetool <- paste0(script_dir,'/../tree/treetools.r')

suppressWarnings({
    source(util)
    source(treetool)
    library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--multigene_table", help="Path to multigene table")
parser$add_argument("--bowtie_aln", help="Path to bowtie alignment file")
parser$add_argument("--metadata", help="Path to metadata")
parser$add_argument("--sf_dir", help="Path to StrainFinder directory")
parser$add_argument("--phylo_enrichment_dir", help="Path to directory containing phylo enrichment results")
parser$add_argument("--output_dir", help="Path to output directory")
parser$add_argument("--genomes", help="List of genomes that we run StrainFinder on")
args <- parser$parse_args()

genomes.use = args$genomes
sf.tips = args$strainfinder_tips
sf.dir = args$sf_dir 

meta = read.metadata.nice(args$metadata)
rownames(meta) = meta$sample

# load bacterial metadata 
bact = ffread(args$multigene_table, as.dt=T)

# gene counts
# -----------
# ab.g = gene counts matrix [5653 samples x 1987 genes]
# note: contains redundant genes
ab.g = ffread(args$bowtie_aln, row.names=1, header=T)
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
                             
# helper function
load_sf_res = function(gi){
    dir.path = sf.dir
    ifn1 = sprintf('%s/%s.sf.long.txt.3',dir.path, gi)
    ifn2 = sprintf('%s/%s.sf_cons.tab',dir.path, gi)
    u = data.frame(fread(ifn1)[,V1:=strip.dx_labs(V1)][],row.names = 1)
    lb = c(fread(ifn2, select=1)[[1]],'other')
    colnames(u) = paste(gi,lb,sep='.')
    u
}

# get genome that we have strain abundance estimates
gis = readLines(genomes.use)

# load phylo enrichment
ph = do.call(rbind, sapply(gis, function(gi) fread(sprintf('%s/%s.phylo.enrich.full.txt', args$phylo_enrichment_dir, gi)) ,simplify=F))
ph = ph[test == 'nice_dxIBD' | is.na(test)]
ph[, padj:=p.adjust(pval, 'fdr'), .(name)]
ph[, id := do.call(paste, c(.SD[,1:2], sep = "."))]
ph[, dir :=ifelse(coef>0,'D','H')]
setkey(ph,id)
                           
use = sapply(gis, function(gi){
        u = load_sf_res(gi)
        ph[colnames(u),on=.(id)][padj<0.05][!is.na(dir)][,lenu(dir)==2]})
                           
gis = intersect(gis[use], colnames(tpm.b))
                           
# load data and scale
v = sapply(gis, function(gi){

    # load data
    u = load_sf_res(gi)
    id2type = ph[colnames(u),on=.(id)][!is.na(dir)]
    
    # align data
    i = intersect(rownames(u), rownames(tpm.b))
    u = u[i,]
    h = id2type[dir == 'H'][which.min(pval)]
    d = id2type[dir == 'D'][which.min(pval)]
    
    # store data
    H.orig = u[, h$id]
    D.orig = u[, d$id]
    H.padj = h$padj
    D.padj = d$padj
    H.id = h$id
    D.id = d$id
    ui = data.frame(h.id=H.id, h.padj=H.padj, H.orig=H.orig, d.id=D.id, d.padj=D.padj, D.orig=D.orig)
    
    # normalize
    ui$h.norm = ui$H.orig/(ui$H.orig + ui$D.orig)
    ui$d.norm = ui$D.orig/(ui$H.orig + ui$D.orig)
    
    # abundance
    ui$H.ab = log2( (ui$H.orig * tpm.b[i, gi]) + 1)
    ui$D.ab = log2( (ui$D.orig * tpm.b[i, gi]) + 1)
    
    # add index
    ui$id = i
    ui = data.table(merge(ui, long, by='id',all.x=T))[!is.na(day)][,name:=gi]
    ui[,fcal0:=fcal]
    ui[,fcal:=log2(fcal+1)]
    ui[complete.cases(ui)]
    
},simplify=F)
v = do.call(rbind,v)
v$cohort = meta[v$id, on=.(id),cohort]       
v$dx = meta[v$id, on=.(id),nice_dx]
                           
# if there are replicates select the one with the largeset fcal 
v = v[,.SD[which.max(fcal)],.(name, cohort, subject, day)]
                           
# ---------------------------------
# select time points for experiment
# ---------------------------------

# filtering criteria 
FREQ_CUTOFF = 0.1
FCAL_CUTOFF = 1.5
TIME_CUTOFF = 14

# helper funciton that goes through all combinations
select_points = function(u){
    
    if(nrow(u) <= 1){return(NULL)}
    
    # fix variables for pairs
    m = c(H.orig='h.fq', D.orig='d.fq', H.sig='h.norm', D.sig='d.norm', H.ab='h.ab', D.ab='d.ab', fcal0='c.raw', fcal='c')
    j = colnames(u) %in% names(m)
    colnames(u)[j] = m[colnames(u)[j]]
    
    # scale measurements
    u$h.zfq = scale(logit(u$h.fq))
    u$d.zfq = scale(logit(u$d.fq))
    u$h.zab = scale(u$h.ab)
    u$d.zab = scale(u$d.ab)
    u$c.z   = scale(u$c)
    
    # enumerate all pairs
    pairs = as.data.table(data.frame(expand.grid(1:nrow(u), 1:nrow(u))))
    colnames(pairs) = c('i0', 'iF')
    
    # get measurements at t=0
    u0 = u[pairs$i0][,.(cohort, dx, day, id, h.fq, d.fq, h.zfq, d.zfq, h.norm, d.norm, h.ab, d.ab, h.zab, d.zab, c.raw, c, c.z, h.id, d.id, h.padj, d.padj)]
    colnames(u0) = gsub('id', 'id0', gsub('day', 't0', gsub('^c$', 'c0', gsub('\\.', '0.', colnames(u0)))))
    
    # get measurements at t=f
    uf = u[pairs$iF][,.(day, id, h.fq, d.fq, h.zfq, d.zfq, h.norm, d.norm, h.ab, d.ab, h.zab, d.zab, c.raw, c, c.z, h.id, d.id, h.padj, d.padj)]
    colnames(uf) = gsub('id', 'idf', gsub('day', 'tf', gsub('^c$', 'cf', gsub('\\.', 'f.', colnames(uf)))))

    # combine all matrices
    pairs = cbind(cbind(pairs, u0), uf)
    
    # select optimal pair
    pairs = pairs[tf > t0]
    pairs[, delta_c := cf - c0]
    pairs[, delta_t := tf - t0]
    pairs = pairs[pmin(h0.fq, d0.fq) > FREQ_CUTOFF][delta_t > TIME_CUTOFF][pmax(c0.raw, cf.raw) > FCAL_CUTOFF*pmin(c0.raw, cf.raw)]
    pairs = pairs[order(-abs(delta_c))]
    pairs[1]
    
}

# actually subset data
vi = v[,select_points(.SD),.(name, subject)]
vi = vi[complete.cases(vi)]
vi[,c("i0","iF") := NULL][, dir:= sign(delta_c)]

# rename and fix variables
x = vi
x = x[,.(name, cohort, subject, dx, id0, idf, t0, tf, delta_t,
         h0.id0, h0.padj, d0.id0, d0.padj,
         h0.ab, d0.ab, hf.ab, df.ab, h0.zab, d0.zab, hf.zab, df.zab,
         h0.fq, d0.fq, hf.fq, df.fq, h0.zfq, d0.zfq, hf.zfq, df.zfq,
         h0.norm, d0.norm, hf.norm, df.norm,
         c0.raw, cf.raw, c0, cf, c0.z, cf.z)]
x = unique(x)

x = setnames(x, old=c('h0.id0', 'h0.padj', 'd0.id0', 'd0.padj'), new=c('h.id', 'h.padj', 'd.id', 'd.padj'))

col.use = c('t0', 'tf',
     'h0.ab', 'd0.ab', 'hf.ab', 'df.ab',
     'h0.zab', 'd0.zab', 'hf.zab', 'df.zab',
     'h0.fq', 'd0.fq', 'hf.fq', 'df.fq',
     'h0.zfq', 'd0.zfq', 'hf.zfq', 'df.zfq',
     'h0.norm', 'd0.norm', 'hf.norm', 'df.norm'
    )

x[,(col.use) := lapply(.SD, as.numeric), .SDcols=col.use]

# calculate deltas
# ----------------
x[, delta_h.fq := hf.fq - h0.fq]
x[, delta_d.fq := df.fq - d0.fq]

x[, delta_h.zfq := hf.zfq - h0.zfq]
x[, delta_d.zfq := df.zfq - d0.zfq]

x[, delta_c := cf - c0]

x[, cdir := ifelse(c0 < cf, 1, -1)]

# write to file
# --------------
fwrite(x, paste0(args$output_dir,'/strain.competition.csv'), row.names=TRUE)