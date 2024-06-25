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
parser$add_argument("--uhgg_dir")
parser$add_argument('--info')
parser$add_argument('--phylo_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

tree_dir = args$phylo_dr

# requires:
#    '/broad/smillie-data/proj/mgxevo/gene_strain.new/stats.summary.txt'
#       - genereated using: summarize_trees.r
#---------
# metadata
#---------
meta = ffread(paste0(args$uhgg_dir,'/genomes-all_metadata.nice.tsv'))
rownames(meta) = meta$Genome

df = ffread(args$info,as.dt=T)
# only look at nodes that >=5 references below them
df = df[ref.below>=5]

# first look at genomes that have a both a sig IBD/HC node
dg = df[padj<0.05, if(all(c(-1,1) %in% sign(coef))){.SD},.(genome)]
# get all combinations of disjoint nodes between IBD/HC
nds = sapply(unique(dg$genome), function(gi){
        nd = dg[genome==gi][coef>0]$new.nodes
        nh = dg[genome==gi][coef<0]$new.nodes
        tree = read.tree.nice(sprintf('%s/%s.ref.tree',tree_dir, gi)))
    
        nd = setNames(phangorn::Descendants(tree, nd), paste0('node',nd))
        nh = setNames(phangorn::Descendants(tree, nh), paste0('node',nh))

        td = sapply(nd, function(vi) tree$tip.label[vi[vi<=Ntip(tree)]],simplify=F)
        th = sapply(nh, function(vi) tree$tip.label[vi[vi<=Ntip(tree)]],simplify=F)

        # rows = health, cols = disease
        disjoint = do.call(cbind, sapply(td, function(a) do.call(rbind, sapply(th, function(b) length(intersect(a,b))==0, simplify=F)),simplify=F))
        rownames(disjoint) = names(nh)
        colnames(disjoint) = names(nd)                                                   
        if(!any(disjoint)){return(NA)}
        ai = arrayInd(which(disjoint),dim(disjoint))
        ai = cbind(rownames(disjoint)[ai[,1]], colnames(disjoint)[ai[,2]])
        colnames(ai) = c('H','D')
        ai = data.frame(gsub('node','',ai))
        ai = cbind(genome=gi, data.frame(sapply(ai, as.integer)))
        ai},simplify=F)
nds = drop.na(nds)
nds = data.table(do.call(rbind,nds))
nds[,id:=1:nrow(nds)]
                                                                 
dh = merge(nds, dg[,.(genome, D=new.nodes, Rd = ref.below, D.pval = padj, D.coef = coef)], by =c('genome','D'),all.x=T)[order(id)]
dh = merge(dh,  dg[,.(genome, H=new.nodes, Rh = ref.below, H.pval = padj, H.coef = coef)], by =c('genome','H'),all.x=T)[order(id)][,id:=NULL][]
dh = dh[,.SD[order(-Rd,D.pval,-Rh,H.pval)],.(genome)]
                                                                 
dg = df[!(genome %in% dh$genome)]
dg = dg[, if(all(c(-1,1) %in% sign(coef))){.SD},.(genome)]
dg = dg[genome %in% dg[,.SD[padj<0.05][which.min(padj)],.(genome)]$genome]

nds = sapply(unique(dg$genome), function(gi){
        nd = dg[genome==gi][coef>0]$new.nodes
        nh = dg[genome==gi][coef<0]$new.nodes
        tree = read.tree.nice(sprintf('%s/%s.ref.tree',tree_dir, gi))

        nd = setNames(phangorn::Descendants(tree, nd), paste0('node',nd))
        nh = setNames(phangorn::Descendants(tree, nh), paste0('node',nh))


        td = sapply(nd, function(vi) tree$tip.label[vi[vi<=Ntip(tree)]],simplify=F)
        th = sapply(nh, function(vi) tree$tip.label[vi[vi<=Ntip(tree)]],simplify=F)

        # rows = health, cols = disease
        disjoint = do.call(cbind, sapply(td, function(a) do.call(rbind, sapply(th, function(b) length(intersect(a,b))==0, simplify=F)),simplify=F))
        rownames(disjoint) = names(nh)
        colnames(disjoint) = names(nd)                                                   
        if(!any(disjoint)){return(NA)}
        ai = arrayInd(which(disjoint),dim(disjoint))
        ai = cbind(rownames(disjoint)[ai[,1]], colnames(disjoint)[ai[,2]])
        colnames(ai) = c('H','D')
        ai = data.frame(gsub('node','',ai))
        ai = cbind(genome=gi, data.frame(sapply(ai, as.integer)))
        ai},simplify=F)
nds = drop.na(nds)
nds = data.table(do.call(rbind,nds))
nds[,id:=1:nrow(nds)]
dk = merge(nds, dg[,.(genome, D=new.nodes, Rd = ref.below, D.pval = padj, D.coef = coef)], by =c('genome','D'),all.x=T)[order(id)]
dk = merge(dk,  dg[,.(genome, H=new.nodes, Rh = ref.below, H.pval = padj, H.coef = coef)], by =c('genome','H'),all.x=T)[order(id)][,id:=NULL][]
dk = dk[,.SD[order(-Rd,D.pval,-Rh,H.pval)],.(genome)][,dx.isIBD:=abs(D.pval-pmin(D.pval,H.pval))<1e-12,.(genome)][]
                                                                 
df.use = rbind(dh[,.SD,.(genome)][,dx.isIBD:=NA][], dk[,.SD,.(genome)])
df.use[,test:=ifelse(is.na(dx.isIBD),'H/D', NA)]
df.use[is.na(test),test:=ifelse(dx.isIBD, 'D/other(H)','H/other(D)')]
df.use[,dx.isIBD:=NULL]
                                                                 
# for each tree get a list of refs that we can look at 
rr = sapply(unique(df.use$genome), function(gi){
    ifn = sprintf('%s/%s.ref.tree',tree_dir, gi)
    tree = read.tree.nice(ifn)
    nodes = unique(c(df.use[genome==gi]$H, df.use[genome==gi]$D))
    nodes = setNames(phangorn::Descendants(tree, nodes), paste0('node',nodes))
    tips = sapply(nodes, function(vi) tree$tip.label[vi[vi<=Ntip(tree)]],simplify=F)
    tips = sapply(tips, function(vi){
        i = grepl('^GUT',vi)
        fix_ref_names(vi[i])})
    tips = sapply(tips, function(a) paste(a,collapse=';'))
    cbind(genome=gi, tips, node = gsub('node','',names(tips)))
},simplify=F)
rr = data.table(do.call(rbind,rr))
rr[,node:=as.integer(node)]
df.use = merge(df.use,rr[,.(genome, D.tips=tips, D=node)], by=c('genome', 'D'),all.x=T)
df.use = merge(df.use,rr[,.(genome, H.tips=tips, H=node)], by=c('genome', 'H'),all.x=T)
df.use[,nd.group:=as.numeric(factor(paste(genome, test, D.tips, H.tips)))]
df.use[,D.tips:=NULL][,H.tips:=NULL]

fwrite(df.use, file = sprintf('%s/ids.all.txt', args$outdir), sep='\t', quote=F, row.names = F)

df.use = df.use[,.SD[order(-Rd, D.pval,-Rh, H.pval)],.(genome)]
df.use  = df.use[,.SD[which.min(D.pval*H.pval)],.(nd.group)]
df.use2 = df.use[,.SD[which.min(D.pval*H.pval)],.(genome)]

ids = data.frame(df.use2)
rownames(ids) = ids$genome
                                                                 
ri = sapply(ids$genome, function(gi){
    # read tree 
    ifn = sprintf('%s/%s.ref.tree',tree_dir, gi)
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
    mp},simplify=F)
has.pan = sapply(ri, function(a){mi = unique(a$mg); ifn = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/genes_presence-absence.tsv',args$uhgg_dir,substr(mi,1,13), mi); all(file.exists(ifn))})
ids = ids[which.names(has.pan),]
df.use2 = data.table(ids)
                  
df.s = ffread(args$info,as.dt=T)
df.use2 = merge(df.use2, df.s[,.(genome, H.old=old.nodes, H =new.nodes)], by=c('genome','H'), all.x=T)
df.use2 = merge(df.use2, df.s[,.(genome, D.old=old.nodes, D =new.nodes)], by=c('genome','D'), all.x=T)
fwrite(df.use2, file = sprintf('%s/ids.use.txt', args$outdir), sep='\t', quote=F, row.names = F)