# input: 
#   gi = genome ID from `ids.use.txt`
# requires:
#   *.locus.mapped.csv | *.uhgp.csv = generated using map2uhgp90.r

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
parser$add_argument('--phylo_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

gi = args$genome
uhgg.dir = args$uhgg_dir
tree.dir = args$phylo_dir
# genome 
# ids.use.txt
# trees with refs
# UHGG catalog dir 
# output dir

#--------------
# load metadata
#--------------
# ids = ffread('/broad/smillie-data/proj/mgxevo/gene_strain.new/ids.use.txt')
ids = ffread(args$gene_strain_ids)
rownames(ids) = ids$genome

meta = ffread(paste0(uhgg.dir,'/genomes-all_metadata.tsv'))
rownames(meta) = meta$Genome

#----------------
# build dataframe
#----------------
# read tree 
ifn = sprintf('%s/%s.ref.tree',tree.dir, gi)
tree = read.tree.nice(ifn)

# get reference genomes
D = fix_ref_names(grep('^GUT', extract.clade(tree, ids[gi,'D'])$tip.label, value=T))
H = fix_ref_names(grep('^GUT', extract.clade(tree, ids[gi,'H'])$tip.label, value=T))

# we shouldn't have any common genomes but still
both = intersect(D,H)
D = setdiff(D, both)
H = setdiff(H, both)

# build a nice dataframe to access mgygs
mp = rbind(cbind(gid = D, type='D'), cbind(gid = H, type='H'))
mp = data.frame(cbind(mp, mg=meta[mp[,'gid'],]$mgyg))
rownames(mp) = mp$gid
mp$completeness = meta[rownames(mp),]$Completeness/100

#----------------------------------------    
# build UHGP-90 <-> KEGG counts matricies
#----------------------------------------
#construct genome-gene list
mg = unique(mp$mg)
mm = sapply(mg, function(mi){
    ifn = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/genes_presence-absence_locus.mapped.csv', uhgg.dir, substr(mi,1,13), mi)
    if(file.exists(ifn)){
        genomes = colnames(fread(ifn, header=T, nrow=0))
        j = sort(match(intersect(mp$gid,genomes), genomes))
        mm = bigreadr::big_fread2(ifn, select = j)
        mm  = sapply(mm, function(a) unique(unlist(strsplit(a,'\t'),use.names=F)),simplify=F)
    } else {
        ifn2 = sprintf('%s/uhgg_catalogue/%s/%s/genome/%s.uhgp.txt', uhgg.dir, substr(mi,1,13),mi, mi)
        mm = readLines(ifn2)
        mm = setNames(list(mm),meta[meta$mgyg==mi,]$Genome)
    }
    mm},simplify=F)
mm = unlist(unname(mm), recursive = F)

genes.all = unique(unlist(mm,use.names=F))

U = do.call(cbind, sapply(names(mm), function(nm) genes.all %in% mm[[nm]], simplify=F))
rownames(U) = genes.all
U = Matrix(U, sparse=T)*1
                          
# UHGP-90 --> KEGG
kmap = fread(sprintf('%s/uhgp-90/uhgp-90.ko',uhgg.dir))
kmap = kmap[id %in% rownames(U)]
K = sparse_agg(U[kmap$id,], kmap$ko)
                          
#-----------------------------------------
# locus x genome --> gene symbol x genome
#-----------------------------------------
mg = unique(mp$mg)
mm = sapply(mg, function(mi){
    # load paths
    path = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/',uhgg.dir, substr(mi,1,13), mi)
    ifn1 = paste0(path,'genes_presence-absence.tsv')
    ifn2 = paste0(path,'genes_presence-absence_locus.csv')
    ifn3 = paste0(path,'pan-genome_eggNOG.tsv')
    ifn4 = paste0(path,'genes_presence-absence_locus.mapped.csv')
    
    if(!file.exists(ifn2)) return(NA)
    
    # get relevant cols and read abs/pres matrix
    j = colnames(fread(ifn1, header=T, nrow=0))
    j = c(1, sort(match(intersect(mp$gid, j), j)))
    A = bigreadr::big_fread2(ifn1, select = j)
    A  = data.frame(A ,row.names=1)
    A  = A[rowSums(A)>0,,drop=F]
    A  = Matrix(data.matrix(A), sparse=T)

    # load the gene names + the "alternative" names called by Roary
    gmap = data.table(bigreadr::big_fread2(ifn2, select = 1:2))
    colnames(gmap) = c('gene','roary')
    gmap[,row:=1:nrow(gmap)]
    gmap[,mgyg:=mi]

    # get the locus <-> KEGG mappings
    km = fread(cmd=paste('cut -f1,6', ifn3))
    colnames(km) = c('locus','kegg')
    km = km[kegg!='']

    # map locus <-> gene-symbol (via row number)
    u2i = unique(get_uhgp_index(unique(km$locus), ifn2))[!is.na(gene)][gene!='']

    # merge: kegg <-> gene symbol
    km = merge(km, u2i[,.(locus=gene, row)], by='locus',all=T)[order(row)]
    gmap = merge(gmap, km, by='row', all=T)
    gmap[gmap=='']=NA

    # select gene-symbols:
    # 1. use the default gene names
    # 2. if a gene is unasigned (i.e., group[0-9]) use the "non-unique" name called by Roary
    # 3. if there's no name called by Roary, use the name called by KEGG
    gmap[,symbol:=gene]
    gmap[grepl('group',symbol) & !is.na(roary),symbol:=roary]
    gmap[grepl('group',symbol) & !is.na(kegg),symbol:=kegg]
    gmap[grepl('group', symbol), symbol:=paste(mi,symbol,sep='_')]

    # tidy the symbols and focus on legit gene symbols
    gmap[!grepl('MGYG',symbol), symbol:=gsub('_(.*)$','',symbol)]
    gmap = gmap[grepl('^[a-z]',symbol)][nchar(symbol)>2]
    gmap[,mgyg:=mi]
    # align and agg
    i = intersect(rownames(A), gmap$gene)
    A = A[i,,drop=F]
    gmap = gmap[gene %in% i]
    A = sparse_agg(A[gmap$gene,,drop=F], gmap$symbol)
    list(A=A, gmap=gmap)
}, simplify=F)
mm = mm[!is.na(mm)]

# intialise matrix dim names
cols = unique(unlist(sapply(mm, function(a) colnames(a[['A']]), simplify=F),use.names=F))
rows = unique(unlist(sapply(mm, function(a) rownames(a[['A']]), simplify=F),use.names=F))

G = Matrix(0, nrow=length(rows), ncol=length(cols), sparse=TRUE)
rownames(G) = rows
colnames(G) = cols
                            
for(mt in mm){
    G = G + pad_mtx(Matrix(mt[['A']], sparse=T), rows = rows, cols = cols)
}
gmap = do.call(rbind, sapply(mm, function(a) a[['gmap']], simplify=F))

# save output
res = list(U=U, G=G, K=K, d=mp, kmap=kmap, gmap=gmap)
saveRDS(res,sprintf('%s/%s.mtx.rds',args$outdir,.gi))
