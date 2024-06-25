script_path <- this.path::this.path()
script_dir <- dirname(script_path)

mtx <- paste0(script_dir,'/../../ut/mtx.r')
plot <- paste0(script_dir,'/../../pl/plot.r')
colors <- paste0(script_dir,'/../../pl/colors.r')
singlecell <- paste0(script_dir,'/../single_cell/singlecell.r')

suppressMessages({
require("ggplot2");
library("wesanderson");
require("ggtree");
library("Hmisc")
library("ggbeeswarm");
require("patchwork");
library("phytools");
library("phangorn");
library("stringi");
library("ape");
library("stringr");
library('dplyr');
library('data.table');
library('stringdist');
library(tidyverse);
library(grid);
library(lme4);
library(lmerTest);
library(glmnet);
source(mtx);
source(plot);
source(colors);
# source(singlecell);
})

simple_regression = function(x, y, signif.fig=2, cor.type='spearman',...){
    simple_scatter(x,y,...) + 
    geom_smooth(aes(x=x,y=y),method='lm', color='black') + 
    ggtitle(sprintf('Spearman rho = %s (p=%s)', signif(cor.test(x,y,method=cor.type,exact=F)$estimate,signif.fig), signif(cor.test(x,y,method=cor.type,exact=F)$p.value,signif.fig)))
}

lab2trait.refs = function(lab){
     trait = gsub(".*_", "", lab)
trait = gsub("IBDU", "CD", trait)
trait[grepl('^GUTREF',trait)] = 'Reference'
trait = factor(trait, levels = c("CD", "HC", "UC","Reference"), ordered = T)
trait = setNames(trait, lab)
return(trait)
    }

# filter matrix for DE 
filter_mtx_strain_de = function(M, mp){
    M = as(M, "lMatrix")
    M[M!=0] = 1
    M = M[rowSums(M)!=0, ]
    i = rowSums(sapply(unique(mp$type),
               function(vi){
                   ci = mp[mp$type==vi,]$gid
                   rowSums(M[,ci,drop=F])/length(ci) > .05})) > 0
    M = M[i,]
    M
}
    
# given a list of genes and locus mtx, return the row in which that said gene symbol can be found
# TODO: return col number too (easy if mtx is the locus mtx, harder if it's the uhgp-90 mapped mtx
get_uhgp_index = function(genes, mtx){
    ofn = tempfile()
    writeLines(genes, ofn)
    cmd.row = sprintf("LC_ALL=C fgrep -n -o -f %s %s", ofn, mtx)
    i = fread(cmd=cmd.row, sep=':',header=F)
    colnames(i) = c('row','gene')
    i[,row:=row-1]
    i}

node.disjoint.matrix = function(nodes,tree){
    nodes = setNames(phangorn::Descendants(tree, nodes), paste0('node',nodes))
    tips = sapply(nodes, function(vi) tree$tip.label[vi[vi<=Ntip(tree)]],simplify=F)
    disjoint = do.call(cbind, sapply(tips, function(a){
                    ri = do.call(rbind, sapply(tips, function(b) length(intersect(a,b))==0, simplify=F))
                    rownames(ri) = names(tips)
                    ri},simplify=F))
    colnames(disjoint) = names(tips)
    return(disjoint)}
                                     
uniq = function(x) unique(x)
drop.na = function(x) x[!is.na(x)]
na.replace = function(x, r){y=x; y[is.na(y)] = r; return(y)}

label_decimals = function(x){gsub("\\.0$", "", as.character(x))}

quick_treeplot = function(gn, tree_dir, enrich_dir){
    gi = gn
    df2 = fread(sprintf(paste0(enrich_dir, '/', gi, 
    '.phylo.enrich.full.txt')))
    df2 = df2[test=='nice_dxIBD' | is.na(test)][order(node)]
    df2[,padj:= p.adjust(pval,'fdr')]
    
    ifn1 = sprintf(sprintf(paste0(tree_dir, '/', gi, '.tree')))
    ifn2 = sprintf(sprintf(paste0(tree_dir, '/', gi, '.tree.ref')))
    
    if(!file.exists(ifn2)){
        warning('no reference genomes integrate')
        tree = read.tree.nice(ifn1)
        node.use = df2[padj < 0.05,]$node
        node.use = c(Ntip(tree)+1, node.use)
    } else {
        old.tree = read.tree.nice(ifn1)
        tree = read.tree.nice(ifn2)
        node.use = df2[padj < 0.05,]$node
        node.use = c(Ntip(tree)+1, node.map(node.use, old.tree, tree))
    }
    
    anc = data.frame(do.call(rbind, sapply(phangorn::Descendants(tree, 1:(Ntip(tree) + Nnode(tree))), function(ci){
        tips = tree$tip.label[ci[ci < (Ntip(tree) + 1)]]
        table(lab2trait(tips))[c("CD", "HC","UC")]},simplify=F)))
    
    traits = setNames(tree$tip.label, tree$tip.label)
    lvls = c()
    i = grepl('_(CD|UC|HC)$',traits)
    traits[i] = gsub(".*_", "", traits[i])
    lvls = c(lvls,c("CD", "HC", "UC"))
    
    if(file.exists(ifn2)){
        i = grepl('^GUT',traits)
        traits[i] = 'Reference'
        lvls = c(lvls,'Reference')
    }
    traits[!(traits %in% lvls)] = "null"
    lvls = c(lvls, "null")
    traits = factor(traits, levels =lvls)
    
    pal2 = lvls
    pal2[match( c("CD", "HC", "UC"),lvls)] = set.colors[c(1,2,3)]
    if(file.exists(ifn2)){
        pal2[pal2=="Reference"] = '#000000'
    }
    pal2[pal2=='null'] = NA
    
    pie.scale = 0.4
    tip.alpha = 0.5
    tip.scale = 0.2
    lwidth = 0.1
    shrink = 1
    
    if(!file.exists(ifn2)){
        shrink = 0.5
        pie.scale = 1
        tip.alpha = 0.5
        tip.scale = 1
    }

    pdf(sprintf(paste0(tree_dir, '/', gi, '.tree.pdf')), height=6, width=6)
    plot.tree(tree,
              traits = traits,
              col = pal2,
              anc.color = pal2,
              anc=anc,
              shrink= shrink,
              pie.scale = pie.scale * 2,
              tip.alpha = tip.alpha,
              tip.scale = tip.scale,
              anc.legend = FALSE,
              legend.loc = 'topleft',
              node.use = node.use,
              title.scale = 2,
              lwidth = 0.1)
    
    dev.off()

}

taxa2phylo = function(df){
        rdp = apply(df, 1, paste0, collapse='|')
        ifn = tempfile()
        ofn = tempfile()
        writeLines(rdp,ifn)
        system(sprintf('python /broad/smillie-data/code/tree/rdp2newick.py -i %s -o %s -f normal',ifn, ofn))
        tree = read.tree(ofn)
        tree$edge.length <- rep(1, nrow(tree$edge))
        tree$tip.label = gsub('(.*)\\|','',tree$tip.label)
        return(tree)
}

fix_ref_names = function(a){
    ret = str_extract(a, 'GUTREFGUTGENOME\\d{6}')
    ret = gsub('GUTREFGUTGENOME', 'GUT_GENOME', ret)
    return(ret)}

                          
                          
nm_wrap = function(x, width = NULL){
    if(is.null(width)){
        width = nchar(x)
    }
    if(width + 1 >= nchar(x)){
        return(x)
    } else {
        paste(c(substr(x,1,width),substr(x,width+1,nchar(x))),collapse="-\n")
    }}

all.equal2 = function(...){
    ret = as.logical(all.equal(...))
    ret[is.na(ret)] = FALSE
    return(ret)
}
lazy_labs = function(data, meta=NULL){
    if(is.null(meta)){
        meta = fread('/broad/smillie-data/db/uhgg/genes/ids.multigene.nice.txt',data.table=F)
        rownames(meta) = meta$genome
    }
    lab = setNames(meta[data$name,]$nice.species, data$name)
    lab = lab[order(log10(p.adjust(data$padj, 'fdr')))]
    lab[duplicated(lab)] = NA
    lab = lab[data$name]
}

node.map = function(nodes, old.tree, new.tree){
    sapply(nodes, function(node) getMRCA(new.tree, extract.clade(old.tree, node)$tip.label))
}
parse_treestats.glm.single = function (gn, trait_in = "IBD") 
{
if(trait_in=="IBD.full"){
    # warning('IBD.full = IBD for single gx as 18 DEC 22')
    trait = "IBD"
} else {
    trait = trait_in
   }
phy = read.tree.nice(paste0("/broad/smillie-data/proj/mgxevo/phylo/single/", 
gn, ".tree.boot_collapse"))
anc = do.call(rbind, sapply(1:(Ntip(phy) + Nnode(phy)), function(i) {
if (i <= Ntip(phy)) {
    table(lab2trait(phy$tip.label[i]))[c("CD", "HC", 
        "UC")]
}
else {
    table(lab2trait(extract.clade(phy, i)$tip.label))[c("CD", 
        "HC", "UC")]
}
}, simplify = FALSE))
anc = setColnames(anc, paste0(colnames(anc), ".node"))
anc = cbind(setColnames(anc[rep(Ntip(phy) + 1, each = nrow(anc)), 
, drop = F], gsub("node", "all", colnames(anc))), anc)
anc = as.data.frame(cbind(node = 1:(Ntip(phy) + Nnode(phy)), 
anc))
df = ffread(paste0('/broad/smillie-data/proj/mgxevo/phylo/trees/all/',gn,'.glm.full.txt'))
df$test[is.na(df$test)] = "NA"
df = setDT(df)[order(df$node)]
pd = data.frame(test = c("IBD", "IBD.full", "age", "sexF", 
    "bmi"), coef = NA, pval = NA, note = NA)
bloc = df[df$test == "NA", ]
bloc = unique(bloc)
bloc = bloc[order(node)]
bloc = cbind(bloc[rep(1:nrow(bloc), each = nrow(pd)), c("name", 
    "node")], pd[rep(c(1:nrow(pd)), nrow(bloc)), ])
df = df[df$test != "NA", ]
df = rbind(df, bloc)
df = df[order(node)]
df = df[df$test != "(Intercept)", ]
df$test = gsub("nice_|nice_dx|nice_dx2", "", df$test)
df = df[df$test == trait, ]
stats = merge(anc, df[, c("test", "node", "coef", "pval")], 
    by = "node", all.x = TRUE)
stats[, 9:10] = sapply(stats[, 9:10], as.numeric)
stats$isTip = stats$node <= Ntip(phy)
stats$padj = NA
stats$test[is.na(stats$test)] = "NA"
for (trait.test in setdiff(unique(stats$test), "NA")) {
    i = stats$test == trait.test
    stats[i, ]$padj = p.adjust(stats[i, ]$pval, "fdr")
}
colnames(stats) = gsub("test", "var", colnames(stats))
colnames(stats) = gsub("coef", "estimate", colnames(stats))
if(trait_in=="IBD.full"){
    stats$var = gsub("IBD", "IBD.full", stats$var)
   }
    return(stats)
}
get.anc = function(phy){
        anc = do.call(rbind, sapply(1:(Ntip(phy)+Nnode(phy)),function(i){
        if(i <= Ntip(phy)){
            table(lab2trait(phy$tip.label[i]))[c("CD", "HC", "UC")]
            # c(0,0,0)
        } else {
            table(lab2trait(extract.clade(phy,i)$tip.label))[c("CD", "HC", "UC")]
        }},simplify=FALSE))
    
    anc = setColnames(anc, paste0(colnames(anc),'.node'))
    anc = cbind(setColnames(anc[rep(Ntip(phy)+1,each=nrow(anc)),,drop=F], gsub('node','all',colnames(anc))), anc)
    anc = as.data.frame(cbind(node=1:(Ntip(phy)+Nnode(phy)), anc))
}


select_nodes2 = function(tree, stats_in, inner.pval=0.05, outer.pval=0.05, trait.test='IBD.full', single.gene=FALSE){
    m = impute.metadata()
    if(single.gene & trait.test=='IBD.full'){
        stats = stats_in[stats_in$var == 'IBD',]
    } else {
                stats = stats_in[stats_in$var == trait.test,]
}
    stats = stats[stats$padj < outer.pval & !is.bad(abs(stats$estimate)),]
    stats = stats[order(-abs(stats$estimate)),]
    stats$abs.effect = abs(stats$estimate)

    nodes.use = c(rep(NA, nrow(stats)))
    enriched = c()

    if(nrow(stats)==0){
        return(c(Ntip(tree)+1))
    }

    if(nrow(stats)>0){
        for (i in 1:nrow(stats)) {

            # stats @ node i
            ni = stats[i, ]
            ci = extract.clade(tree, ni$node)$tip.label
            si = length(ci)

            #test each node against its descendants
            di = intersect(getDescendants(tree, ni$node), stats$node)
            sig.child = any(sapply(di, function(j){
                # node j stats
                nj = stats[stats$node == j, ]
                if (ni$abs.effect > nj$abs.effect) {
                    return(FALSE)
                }
                # enrichment test
                na = ni$node
                nb = nj$node
                u = enrichment.between(na, nb, phy = tree, metadata = m, type = 'all')
                u = u[u$var==trait.test,]
                u.pval = ifelse(is.na(u$pval), 1, u$pval)
                if (u.pval < inner.pval) {
                    return(TRUE)
                }
                return(FALSE)
            }))

            # iteratively select most enriched nodes
            if (sig.child) {
                nodes.use[i] = FALSE
                next
            }

            # check other nodes
            for (j in which(nodes.use)) {
                if (i == j) {
                    next
                }
                # stats @ node j
                nj = stats[j, ]
                cj = extract.clade(tree, nj$node)$tip.label
                sj = length(cj)

                # if there's an overlap
                if (length(intersect(ci, cj))) {

                    # enrichment test (ni vs. nj)
                    u = enrichment.between(ni$node, nj$node, phy = tree, metadata = m, type = 'all')
                    u = u[u$var==trait.test,]
                    u.pval = ifelse(is.na(u$pval), 1, u$pval)
                    if (u.pval >= inner.pval) {
                      # if not significant, keep largest node
                      if (si > sj) {
                        nodes.use[i] = TRUE
                        nodes.use[j] = FALSE
                      } else {
                        nodes.use[i] = FALSE
                        nodes.use[j] = TRUE
                      }
                    } else {
                      # if significant, keep most enriched node
                      if (ni$abs.effect > nj$abs.effect) {
                        nodes.use[i] = TRUE
                        nodes.use[j] = FALSE
                      } else {
                        nodes.use[i] = FALSE
                        nodes.use[j] = TRUE
                      }
                    }
                }
            }

            # otherwise, keep node
            if (is.na(nodes.use[i])) {
                nodes.use[i] = TRUE
            }
         }
         enriched = c(enriched, stats[nodes.use,]$node)
    }
    new.nodes = c()
    # then for those nodes that arent enriched iterate over the largest disjoint clades that have at least 10 tips
    df = stats_in[stats_in$var == trait.test & !stats_in$isTip,]
    df$size = sapply(df$node, function(a) Ntip(extract.clade(tree,a)))
    df = df[df$size>10,]
    df = df[order(-df$size),]
    while(nrow(df)>0){
        cu = unique(unlist(sapply(c(new.nodes,enriched), function(a) extract.clade(tree,a)$tip.label)))
        overlap = sapply(df$node, function(ni) length(intersect(extract.clade(tree,ni)$tip.label, cu)))
        df = df[!overlap,]
        if(!sum(!overlap)){
            break
        }
        df = df[order(-df$size),]
        new.nodes = c(new.nodes, df[1,]$node)
    }
    return(c(new.nodes,enriched))
}

scientific_10 = function(x){
    text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))
    text[x==0] = "0"
    text[x==1] = "1"
    text=gsub("1 %\\*% 10", "10", text)
    text=gsub("10\\^01","10",text)
    return(parse(text=text))
}

hush=function(code){
  tmp = suppressMessages({
  # sink(file=NULL,type = "message") # use /dev/null in UNIX
      code
  })
  return(tmp)
}

strip.dx_labs = function(a) gsub('_(UC|CD|HC|IBDU)$','',a)

plot_grid2 = function(..., title=NULL){
    plot_row = plot_grid(...)
    # now add the title
    title.draw = ggdraw() + 
        draw_label(title, fontface = 'bold', x = 0, hjust = 0) +
        theme(plot.margin = margin(0, 0, 0, 7))
    plot_grid(title.draw, plot_row, ncol = 1, rel_heights = c(0.1,1))
}

flatten = function(X){
    unname(unlist(c(X)))
}

setRownames = function(obj, nm){
    rownames(obj) = nm
    return(obj)
}

setColnames = function(obj, nm){
    colnames(obj) = nm
    return(obj)
}

which.names = function(X){
    names(which(X))
}
# re-root tree, B, to align with a rooted tree, A
#------------------------------------------------
align.root = function(A, B){
    tree1 = A
    tree2 = B
    
    # append node labels to tree2 
    f = fortify(tree2)
    f[!f$isTip,]$label = paste0('node_',1:Nnode(tree2))
    tree2 = as.phylo(f)
    tree2.init = tree2
    tree2 = drop.tip(tree2, setdiff(tree2$tip.label, tree1$tip.label))

    # triplet distance 
    ## capture output to hide reroot printing a message
    d = sapply(1:nrow(tree2$edge), function(nd) triplet.dist(tree1,reroot(tree2,nd)))
    i = which.min(d)

    # find the root node for the _original_ tree
    nd = which(fortify(tree2.init)$label == fortify(tree2)[i,]$label)
    res = reroot(B, nd)
    
    # align tree tips
    tree1 = read.tree(text=write.tree(ladderize(tree1)))
    res   = read.tree(text=write.tree(ladderize(res)))
    return(list(tree1, res))
}
               
triplet.dist = function (A, B) 
{
    fn = paste0(stringi::stri_rand_strings(1, 64), c("_A", "_B"), 
        ".nh")
    write.tree(phy = A, fn[1])
    write.tree(phy = B, fn[2])
    cmd = paste("/broad/smillie-data/home/akumbhar/bin/triplet_dist", fn[1], fn[2], "2>&1; rm ",fn[1], fn[2])
    res = system(cmd, intern = TRUE)
    if(any(is.na(strtoi(res)))){
        # print the error message from triplet_dist
        stop(paste0(res,collapse='\n'))
    }
    return(strtoi(res))
}

quart.dist = function (A, B) 
{
    fn = paste0(stringi::stri_rand_strings(1, 64), c("_A", "_B"), 
        ".nh")
    write.tree(phy = A, fn[1])
    write.tree(phy = B, fn[2])
    cmd = paste("triplet_dist", fn[1], fn[2], "2>&1; rm ",fn[1], fn[2])
    res = system(cmd, intern = TRUE)
    if(any(is.na(strtoi(res)))){
        # print the error message from triplet_dist
        stop(paste0(res,collapse='\n'))
    }
    return(strtoi(res))
}

# are a given set of nodes disjoint
node.overlap = function(nodes, tree){
    m = crossprod(table(do.call(rbind, sapply(nodes, 
           function(a){data.frame(lab=extract.clade(tree, a)$tip.label, strain = a)}, simplify=FALSE))))
    !all(m[!diag(nrow(m))] == 0)
}

node.disjoint = function(nodes, tree){
    !node.overlap(nodes,tree)
}

table2file = function(x, file, sep='\t'){
    write.table(x,file = file,row.names = F,sep=sep,quote=F) 
}

len = function(x){
    return(length(x))
}

lenu = function(x){
    return(length(unique(x)))
}
symdiff <- function( x, y) { setdiff( union(x, y), intersect(x, y))}

read.tree.nice = function(file, reroot=FALSE){
    phy = read.tree(file)
    if(reroot){
        phy = midpoint.root(phy)
    }
    phy$tip.label = gsub('Indeterminate$', 'IBDU', phy$tip.label)
    phy$tip.label = gsub("IBDU", "CD", phy$tip.label)
    phy$tip.label = gsub("Healthy", "HC", phy$tip.label)
    return(phy)
}

# OUTPUT nodes whose tips DO NOT vary by at least alpha
node.similarity = function(tree, d, alpha = 99.5) {
  sapply(node, function(a) {
      T = extract.clade(tree, a)$tip.label
      return(max(d[T, T]) < 1 - (alpha/100))
  }, simplify = T)
}

labnode = function(x,color='cyan',opacity=0.5,cex=0.25){
    nodelabels(node = x, 
               pie = matrix(rep(c(0,1),length(x)), nrow=length(x),byrow=T), 
               piecol=alpha(color,opacity), 
               cex = cex)
    }

labtip = function(x,color='cyan',opacity=0.5,cex=0.25){
    tiplabels(tip = x, 
               pie = matrix(rep(c(0,1),length(x)), nrow=length(x),byrow=T), 
               piecol=alpha(color,opacity), 
               cex = cex)
    }

get.direct_child = function(x, node){
    v = x$edge[x$edge[, 1] == node, 2]
    return(x$tip.label[v[v<=Ntip(x)]])
}
get.directChild = function(tree,node){
    get.direct_child(tree, node)
}

get.terminalNodes = function(tr){
    ix = sapply(as.integer(gsub('@.*','',tr$node.label)), function(a){
        if(length(extract.clade(tr,a)$node.label)==1){
            return(a)
        } else{
            return(NA)
        }
    },simplify=T)
    ix = ix[!is.na(ix)]
    return(ix)
}


# is ni a parent of the nodes listed in nj
is.parent = function(ni,nj, tree){
    children = as.integer(names(extract.clade(tree, ni)$node.label))
    x=sapply(nj, function(a){a %in% children})
    if(!length(x)){
        return(FALSE)
    } else{
        return(x)
    }
}
# is ni a child of the nodes listed in nj
is.child = function(ni,nj, tree){
   x=sapply(nj, function(a){ni %in% as.integer(names(extract.clade(tree, a)$node.label))})
    if(!length(x)){
            return(FALSE)
        } else{
            return(x)
        }
}

get.children = function(node, tree){
    return(as.integer(names(extract.clade(tree,node)$node.label)))
}

tab2trait = function(cluster, tree, numeric=TRUE){
    cluster = data.frame(cluster)
    rownames(cluster) = cluster$tip.label
    g = setNames(cluster[,!grepl('tip.label',colnames(cluster))],cluster$tip.label)
    # g = factor(as.character(g), levels=str_sort(unique(g),numeric=numeric), ordered=T)
    g.bad = names(g[is.bad(g)])
    
    x = setNames(c(rep('NA',Ntip(tree))), tree$tip.label)
    ix = intersect(names(g),names(x))
    x[ix] = g[ix]
    x[g.bad] = 'NA'
    x = factor(x, levels=str_sort(unique(x),numeric=numeric), ordered=T)
    x = x[tree$tip.label]
    return(x)
}

node2tab = function(nodes, tree){
    df = data.frame(tip.label = tree$tip.label, cluster='NA')
    rownames(df) = df$tip.label
    df = df[tree$tip.label,]
    for(i in 1:length(nodes)){
        df[extract.clade(tree,nodes[i])$tip.label,'cluster'] = i
    }
    
    return(df)
}

node2trait = function(nodes, tree){
#     df = data.frame(tip.label = tree$tip.label, cluster='NA')
#     rownames(df) = df$tip.label
#     df = df[tree$tip.label,]
#     for(i in 1:length(nodes)){
#         df[extract.clade(tree,nodes[i])$tip.label,'cluster'] = i
#     }
    df = node2tab(nodes,tree)
    return(tab2trait(df, tree))
}

get.dist = function(fasta, metric = "hamming", ...) {
  seq = read.dna(fasta, format = "fasta")
  if (metric == "hamming") {
    d = dist.hamming(seq, ...)
  } else if (metric == "p") {
    d = dist.p(seq, ...)
  } else {
    d = dist.ml(seq, model = metric, ...)
  }
  d = as.data.frame.matrix(as.matrix(d))
  rownames(d) = gsub("IBDU", "CD", rownames(d))
  colnames(d) = gsub("IBDU", "CD", colnames(d))
  return(d)
}

lab2trait = function(lab) {
  trait = gsub(".*_", "", lab)
  trait = gsub("IBDU", "CD", trait)
  extra = setdiff(unique(trait),c("CD", "HC", "UC"))
  trait = factor(trait,
                 levels = c(c("CD", "HC", "UC"),extra),
                 ordered = T)
  trait = setNames(trait, lab)
  return(trait)
}

# based on `as.polytomy` from ggtree note that if using node numbers they need
# to be offset by Ntip(tree)
collapse.nodes = function(tree, nodes) {
  feat = 1:tree$Nnode + Ntip(tree)
  if (is.logical(nodes)) {
    idx = which(nodes)
    edge_idx = match(feat[idx], tree$edge[, 2])
  } else if (is.numeric(nodes)) {
    edge_idx = match(nodes, tree$edge[, 2])
  }
  tree$edge.length[edge_idx] = 0
  res = di2multi(tree)
  
  newlab = do.call(rbind, str_split(res$node.label, "@"))
  newlab[, 1] = node.idx(res)
  newlab = apply(newlab, 1, paste, collapse = '@')
  res$node.label = newlab
  return(res)
}

tab2node = function(df, tree){
    sapply(unique(df$cluster), function(a){setNames(getMRCA(phy = tree, tip= df[df$cluster == a,]$tip.label),a)})
}    
plot.tree = function(tree,
                     traits = NULL,
                     stats = NULL,
                     title = "",
                     title.scale = 1,
                     trait.title = "Legend",
                     anc = NULL,
                     anc.color = set.colors[c(1, 3, 2)],
                     anc.title = 'Diagnosis',
                     anc.legend = TRUE,
                     anc.legend.loc = 'topleft',
                     gen.scale = NULL,
                     pie.scale = 1,
                     tip.scale = 1,
                     tip.offset = 1,
                     node.use = NULL,
                     adj.pvals = NULL,
                     pval.correction = 'fdr',
                     min.pval = 0.01,
                     tip.labels = FALSE,
                     tip.traits = TRUE,
                     color = set.colors[c(1, 3, 2)],
                     trait.legend = TRUE,
                     legend.loc = "topright",
                     legend.yspace = 1,
                     anc.legend.inset = 0.05,
                     legend.xspace =1,
                     legend.landscape = FALSE,
                     tip.alpha = 0.8,
                     pie.alpha = 1,
                     na.color = '#cccccc',
                     align.tip.label= FALSE,
                     lty = 1,
                     trait.legend.inset = 0.05,
                     lazy = FALSE,
                     vlazy = FALSE,
                     phylo.type = "fan",
                     lwidth=1, 
                     shrink = 1,
                     font.style = 1,
                     font.scale = 1,
                     ring.node = NULL,
                     ring.color = set.colors[-c(1, 3, 2)], 
                     ring.offset = 0.05, 
                     ring.thickness = 10,
                     ring.legend = TRUE,
                     ring.inset = 0,
                     ring.title = '',
                     ring.legend.loc = 'bottomright',
                     ladderize.tree = FALSE) {
  phylo = tree
  if(ladderize.tree){
      phylo = ladderize(phylo)
  }
  stats = add.tipData(stats, phylo)
  if(is.null(gen.scale)){dx = min(1,max(nodeHeights(phylo)))} else {dx = gen.scale}
  if (lazy & is.null(stats)) {
    stats = calc.treestats2(phylo)
  }
  if (lazy) {
    if(is.null(traits)){
        traits = lab2trait(phylo$tip.label)
    }
    if (!is.null(stats)) {
      anc = stats2anc(stats)
      adj.pvals = stats2pval(stats,method=pval.correction)
    }
    trait.title = "Diagnosis"
  }
  if (!is.null(traits)){
      if(length(levels(traits))>3 & length(color)==3){
          if(all(color == set.colors[c(1, 3, 2)])){
          color = set.colors
        }
      }
    # any bad or emmpty values
    if(any(grepl(pattern = 'na|null|nan',x=levels(traits),ignore.case=TRUE))){
        color=append(color,na.color, after=which(grepl(pattern = 'na|null|nan',x=levels(traits),ignore.case=TRUE))-1)
    }
        col = color[traits]

  }
  if (!is.null(stats) & is.null(adj.pvals)) {
      adj.pvals = stats2pval(stats,method=pval.correction)
    }
  if (!is.null(anc)) {
    if(is.null(node.use)){
        node.use = 1:phylo$Nnode + Ntip(phylo)
        if (!is.null(adj.pvals))
          node.use = c(length(phylo$tip.label) + 1, which(adj.pvals <= min.pval))
    }
  }
  if(!as.logical(length(tip.labels)-1) & !as.logical(tip.labels[1])){ phylo$tip.label[] = ""}
  # if (!as.logical(tip.labels)){
  #     phylo$tip.label[] = ""
  # }
  if(any(is.na(as.logical(tip.labels))) & (length(tip.labels) == length(phylo$tip.label))){
      phylo$tip.label = tip.labels
  }
  if (!is.null(traits)) {
    plot(
      phylo,
      type = phylo.type,
      tip.color = col,
      cex = font.scale,
      cex.main = title.scale,
      main = title,
      edge.width=lwidth,
      align.tip.label = align.tip.label,
      lty= lty,
      font = font.style,
      x.lim = c(-shrink*1.5*max(nodeHeights(tree)), shrink*1.5*max(nodeHeights(tree)))
    )
    if(!is.null(ring.node)){
        cladelab(ring.node, phylo, col=ring.color, offset = ring.offset, thickness = ring.thickness)
    }
  } else {
    plot(phylo,
         font = font.style,
         type = phylo.type,
         cex = font.scale,
         cex.main = title.scale,
         main = title,
         edge.width=lwidth,
         align.tip.label = align.tip.label,
         lty= lty,
          x.lim = c(-shrink*1.5*max(nodeHeights(tree)), shrink*1.5*max(nodeHeights(tree)))
        )
        if(!is.null(ring.node)){
            cladelab(ring.node, phylo, col=ring.color, offset = ring.offset, thickness = ring.thickness)
        }
  }
  if (trait.legend & !is.null(traits)) {
       par(xpd=TRUE)
    legend(
      legend.loc,
      inset = trait.legend.inset,
      title = trait.title,
       bty = "n",
      y.intersp = legend.yspace,
      x.intersp = legend.xspace,
      horiz= legend.landscape,
      legend = levels(traits)[!is.bad(levels(traits))],
      text.width = c(rep(0, length(levels(traits)[!is.bad(levels(traits))]))),
      fill = color[1:length(levels(traits))][!is.bad(levels(traits))]
    )
  }
  if (tip.traits & !is.null(traits)) {
    tip.labs = t(as.data.frame.matrix(table(
      data.frame(trait = traits, tips = names(traits))
    )))

    par(fg = "transparent")
    tip.labs = tip.labs[tree$tip.label,]
    tiplabels(
      pie = tip.labs[names(traits),],
      piecol = alpha(color[1:length(levels(traits))],
                     tip.alpha),
      cex = tip.scale*dx / 3,
      offset = tip.offset * dx * 0.01
    )
  }
  if (!is.null(anc)) {
    par(fg = "black")
    nodelabels(
      node = node.use,
      pie = anc[node.use,,drop=F],
      cex = pie.scale * dx,
      piecol = alpha(anc.color[1:dim(anc)[2]], pie.alpha)
    )
    if(anc.legend){
    par(xpd=TRUE)
    legend(
      anc.legend.loc,
      inset = anc.legend.inset,
      title = anc.title,
      bty = "n",
      legend = colnames(anc),
      fill = anc.color[1:dim(anc)[2]]
    )}
  }
  if (ring.legend & !is.null(ring.node)) {
      par(xpd=TRUE)
    legend(
      ring.legend.loc,
      inset = ring.inset,
      title = ring.title,
      bty = "n",
      y.intersp = legend.yspace,
      x.intersp = legend.xspace,
      horiz= legend.landscape,
      legend = names(ring.node)[!is.bad(ring.node)],
      text.width = c(rep(0, length(names(ring.node)[!is.bad(names(ring.node))]))),
      fill = ring.color[1:length(names(ring.node))][!is.bad(names(ring.node))]
    )
  }
}

stats2anc = function(stats, norm = FALSE) {
  nix = c("CD", "HC", "UC")
  anc = stats[, paste0(nix,'.node')]
  colnames(anc) = gsub(".node", "", colnames(anc))
  anc = data.frame(anc)
  if (norm) {
    anc = anc / rowSums(anc)
  }
  return(anc)
}

stats2pval = function(stats, method = "fdr") {
  pvals = stats$pval
  if(!('isTip' %in% colnames(stats))){
    stop('No tip data found. Try using `stats = add.tipData(stats, tree)`')
  }
  ip = !stats$isTip
  pvals[ip] = p.adjust(pvals[ip], method)
  return(pvals)
}
stats2odds = function(stats, max.only =FALSE) {
    nix = c("CD", "HC", "UC")
    odds = do.call(rbind, sapply(1:nrow(stats), function(ix) {
        p = data.matrix(stats[ix, paste0(nix,'.node')])
        q = data.matrix(stats[ix, paste0(nix,'.all')])
        u = rbind(p, q - p)
        if (any(is.na(u))) {
            x=data.frame(matrix(rep(NA,length(nix)),nrow=1))
            colnames(x) = sapply(1:length(nix), function(a) paste(nix[-a],collapse='vs'))
            return(x)
        } else {
            x = sapply(1:length(nix), function(a) {
                v = u[,-a]
                colnames(v) = gsub(".all|.node", "", colnames(v))
                x = fisher.test(v)
                df = x$estimate
                names(df) = paste(colnames(v), collapse = "vs")
                return(df)
            }, simplify = T)
            return(x)
        }
    }, simplify = FALSE))
    colnames(odds) = paste0("odds.", colnames(odds))  
    max.odds = apply(odds,1,function(a){
        x = abs(log(a))
        x = x[!is.bad(x)]
        if(length(x)<1){
            return(NA)
        } else {
            return(max(x))
        }
    })
    
    if(max.only){
        return(max.odds)
    }
    return(cbind(odds,'abs.log.odds'=max.odds))
}


node.idx = function(tree) {
  return(1:tree$Nnode + Ntip(tree))
}

attach.nodestats = function(tree, stats=NULL) {
  phylo = remove.nodestats(tree)
#   if (all(grepl('@', phylo$node.label))) {
#     warning('Existing data found. Overwriting.')
#     phylo = remove.nodestats(phylo)
#   }
  if(is.null(stats)){
    stats=calc.treestats(phylo)
  }
  phylo$node.label = apply(stats[!stats$isTip,],1,function(x) paste(x,collapse='@'))
  return(phylo)
}

remove.nodestats = function(tree) {
  if (all(!grepl('@', tree$node.label))) {
    return(tree)
  } else {
    tree$node.label = get.nodestat(tree)$label
    return(tree)
  }
}
                           
get.nodestat = function(tree, nix = NULL) {
  if (length(tree$node.label) < 1) {
    stop('No label found')
  }
  if (all(!grepl('@', tree$node.label))) {
    stop('No attached data found. Did you forget to `attach.nodestats`?')
  }
  cols=c("node", "isTip","label", "pval", "padj", "CD.node", "HC.node", "UC.node", "odds.HCvsUC", "odds.CDvsUC", "odds.CDvsHC", "odds.max")
  res = fread(paste(tree$node.label,collapse = '\n'),sep = '@', header = F)
  colnames(res) = cols
  
  if(!is.null(nix)){
    res= res[res$node %in% nix, ]
  } 
  return(res)
}


#### GET USEFUL SPECIES METADATA

                 
is.bad = function(x){
    if(!length(x)) return(TRUE)
    ret = sapply(x, function(x){
        return(is.na(x) | is.infinite(x) | is.null(x) | !length(x) | x=='' | tolower(x) == 'na' | tolower(x) == 'nan' | tolower(x) == 'null')
    })
    return(ret)
}
                 
get.nodesBelow = function(nj,tr){
    return(as.integer(names(extract.clade(tr, nj)$node.label)))
}
                      
select_nodes = function(tree, stats_in, inner.pval = 0.05, outer.pval = 0.05) {
    stats = stats_in[stats_in$padj < outer.pval & !is.bad(stats_in$abs.log.odds),]
    stats = stats[order(-stats$abs.log.odds),]
    nodes.use = c(rep(NA, nrow(stats)))
    nix = c("CD.node", "HC.node", "UC.node")
    
    enriched = c()
    if(nrow(stats)==0){
        return(c(Ntip(tree)+1))
    }
    if(nrow(stats)>0){
        for (i in 1:nrow(stats)) {

            # stats @ node i
            ni = stats[i, ]
            ci = extract.clade(tree, ni$node)$tip.label
            si = length(ci)
        
            #test each node against its descendants
            di = intersect(getDescendants(tree, ni$node), stats$node)
            sig.child = any(sapply(di, function(j){
                # node j stats
                nj = stats[stats$node == j, ]
                if (ni$abs.log.odds > nj$abs.log.odds) {
                    return(FALSE)
                }
                # fisher test
                na = ni[, nix]
                nb = nj[, nix]
                u = fisher.test(rbind(nb,na-nb))
                if (u$p.value < inner.pval) {
                    return(TRUE)
                }
                return(FALSE)
            }))
        
            # iteratively select most enriched nodes
            if (sig.child) {
                nodes.use[i] = FALSE
                next
            }
        
            # check other nodes
            for (j in which(nodes.use)) {
                if (i == j) {
                    next
                }
                # stats @ node j
                nj = stats[j, ]
                cj = extract.clade(tree, nj$node)$tip.label
                sj = length(cj)

                # if there's an overlap
                if (length(intersect(ci, cj))) {

                    # fisher test (ni vs. nj)
                    if (si > sj) {
                      na = ni[, nix]
                      nb = nj[, nix]
                    } else {
                      na = nj[, nix]
                      nb = ni[, nix]
                    }
                    u = fisher.test(rbind(nb, na - nb))
                    if (u$p.value >= inner.pval) {
                      # if not significant, keep largest node
                      if (si > sj) {
                        nodes.use[i] = TRUE
                        nodes.use[j] = FALSE
                      } else {
                        nodes.use[i] = FALSE
                        nodes.use[j] = TRUE
                      }
                    } else {
                      # if significant, keep most enriched node
                      if (ni$abs.log.odds > nj$abs.log.odds) {
                        nodes.use[i] = TRUE
                        nodes.use[j] = FALSE
                      } else {
                        nodes.use[i] = FALSE
                        nodes.use[j] = TRUE
                      }
                    }
                }
            }

            # otherwise, keep node
            if (is.na(nodes.use[i])) {
                nodes.use[i] = TRUE
            }
        }
        enriched = c(enriched, stats[nodes.use,]$node)
    }
    new.nodes = c()
    # then for those nodes that arent enriched iterate over the largest disjoint clades that have at least 10 tips
    df = stats_in[!stats_in$isTip,]
    df$size = sapply(df$node, function(a) Ntip(extract.clade(tree,a)))
    df = df[df$size>10,]

    while(nrow(df)>0){
        cu = unique(unlist(sapply(c(new.nodes,enriched), function(a) extract.clade(tree,a)$tip.label)))
        overlap = sapply(df$node, function(ni) length(intersect(extract.clade(tree,ni)$tip.label, cu)))
        df = df[!overlap,]
        if(!sum(!overlap)){
            break
        }
        df = df[order(-df$size),]
        new.nodes = c(new.nodes, df[1,]$node)
    }
    
    return(c(new.nodes,enriched))
}

add.tipData = function(df, tree){
    if(is.null(df)){
        return(NULL)
    }
    if('isTip' %in% colnames(df)){
        return(df)
    } else{
        return(data.frame(df, isTip = tree$edge[,2]<Ntip(tree)+1))
    }
}
                         
                         
cladelab=function(nodes, tree, col=set.colors, offset = 0.05, thickness = 10){
    # library(plotrix)
    obj = get("last_plot.phylo", envir = .PlotPhyloEnv)
    h = max(sqrt(obj$xx^2 + obj$yy^2))
    lwd = 10
    lend = 1
    key = unique(nodes)
    for(ix in 1:length(key)){
        d = getDescendants(tree, nodes[ix])
        d = sort(d[d <= Ntip(tree)])
        deg = atan(obj$yy[d]/obj$xx[d]) * 180/pi
        ii <- intersect(which(obj$yy[d] >= 0), which(obj$xx[d] < 0))
        deg[ii] <- 180 + deg[ii]
        ii <- intersect(which(obj$yy[d] < 0), which(obj$xx[d] < 0))
        deg[ii] <- 180 + deg[ii]
        ii <- intersect(which(obj$yy[d] < 0), which(obj$xx[d] >= 0))
        deg[ii] <- 360 + deg[ii]
        plotrix::draw.arc(x = 0, y = 0, radius = (1+ offset) * h, deg1 = min(deg), deg2 = max(deg), lwd = thickness, col = col[ix], lend = 1)
    }
}
                         
lab2trait2 = function (tree.tips, trait = "nice_dx", meta=ffread('/broad/smillie-data/data/mgx/All-IBD/mgx4.csv',row.names = T)){
    lab = gsub("_CD|_UC|_HC|_IBDU", "", tree.tips)
    if (trait == "nice_dx") {
        ret = meta[lab, trait]
        ret = gsub("IBDU", "CD", ret)
        ret = factor(ret, levels = c("CD", "HC", "UC"), ordered = T)
        ret = setNames(ret, tree.tips)
    } else if (trait == 'sex') {
        ret = meta[lab, trait]
        ret = factor(ret, levels = c("M", "F"), ordered = T)
        ret = setNames(ret, tree.tips)
    } else if (trait == 'age_grp') {
        ret = meta[lab, trait]
        ret = factor(ret, levels = c("C", "T", "A", "S"), ordered = T)
        ret = setNames(ret, tree.tips)
    } else {
        ret = factor(meta[lab, trait], ordered = T)
        # test if boolean
        if (length(setdiff(levels(ret), c(TRUE, FALSE))) == 0) {
            ret = factor(meta[lab, trait], levels = c(TRUE, FALSE), 
                ordered = T)
        }
        ret = setNames(ret, tree.tips)
    }
}

stats2odds2 = function(stats, max.only =FALSE, trait.by='nice_dx') {
    nix = setdiff(gsub(".all|.node", "", grep(trait.by, colnames(stats), value = T)), paste0(trait.by,'.pval'))
    nix.boolean = length(setdiff(gsub(paste0(trait.by, "."), "", nix), c(TRUE, FALSE))) == 0

    odds = do.call(rbind, sapply(1:nrow(stats), function(ix) {
        p = data.matrix(stats[ix, paste0(nix, ".node")])
        q = data.matrix(stats[ix, paste0(nix, ".all")])
        u = rbind(p, q - p)
        num.combo = ncol(combn(nix, 2))
        if (any(is.na(u))) {
            x = data.frame(matrix(rep(NA, num.combo), nrow = 1))
            if (num.combo == 1) {
                colnames(x) = paste0(trait.by,'.',gsub(paste0(trait.by,'.'), "",paste(nix, collapse = "vs")))
            }
            else {
                colnames(x) = sapply(1:num.combo, function(a) paste0(trait.by,'.',gsub(paste0(trait.by,'.'), "",paste(combn(nix,2)[,a], collapse = "vs"))))
            }
            return(x)
        }
        else {
            if (num.combo == 1) {
                v = u
                colnames(v) = gsub(".all|.node", "", colnames(v))
                x = fisher.test(v)
                df = x$estimate
                names(df) = gsub(paste0(trait.by,'.'), "",paste(colnames(v), collapse = "vs"))
                names(df) = paste0(trait.by,'.',names(df))
                x = df
            }
            else {
                x = sapply(1:num.combo, function(a) {
                    j = combn(nix,2)[,a]
                    v = u
                    colnames(v) = gsub(".all|.node", "", colnames(v))
                    v = v[,j]
                    x = fisher.test(v)
                    df = x$estimate
                    names(df) = gsub(paste0(trait.by,'.'), "",paste(colnames(v), collapse = "vs"))
                    names(df) = paste0(trait.by,'.',names(df))
                    return(df)
                }, simplify = T)}
        return(x)}
    }, simplify = FALSE))
    colnames(odds) = paste0("odds.", colnames(odds))
    colnames(odds) = gsub(paste0(trait.by,"."),"", colnames(odds))
    max.odds = apply(odds,1,function(a){
        x = abs(log(a))
        x = x[!is.bad(x)]
        ifelse(length(x)<1, NA, max(x))})
    
    if (max.only) {
        return(max.odds)
    }
    ret = cbind(odds, abs.log.odds = max.odds)
    # colnames(ret) = gsub('^nice_dx.', '', colnames(ret))   
    return(ret)
}
     
impute.metadata = function(m = ffread('/broad/smillie-data/data/mgx/All-IBD/mgx4.csv')){
    m$nice_dx = factor(gsub("IBDU", "CD",m$nice_dx), levels=c('HC', 'CD', 'UC'))
    m$nice_dx2 = m$nice_dx
    m$nice_dx = factor(gsub("CD|UC", "IBD",m$nice_dx), levels=c('HC', 'IBD'))
    m$age = as.numeric(m$age)
    m$sex = factor(m$sex, levels=c('M', 'F'))

    if(sum(!is.na(m$age)) < 10){m$age = rnorm(mean=0, sd=1, nrow(m))}
    if(sum(!is.na(m$bmi)) < 10){m$bmi = rnorm(mean=0, sd=1, nrow(m))}
    if(sum(!is.na(m$sex)) < 10){m$sex = sample(c('M', 'F'), nrow(m), replace=T)}
    m$nice_age = scale(ifelse(is.na(m$age), mean(m$age, na.rm=T), m$age))
    m$nice_bmi = scale(ifelse(is.na(m$bmi), mean(m$bmi, na.rm=T), m$bmi))
    m$nice_sex = m$sex
    i = is.na(m$nice_sex)
    m$nice_sex[i] = sample(m$nice_sex[!i], sum(i), replace=T)
    m$nice_sex = factor(m$nice_sex, levels=c('M', 'F'))
    m
}
                                 
                                 
                                 
enrichment.between = function(node1, node2, phy, metadata = impute.metadata(), type='all'){
    tr1 = extract.clade(phy,node1)
    tr2 = extract.clade(phy,node2)

    if(length(intersect(tr1$tip.label, tr2$tip.label))==0){
        stop('Nodes are not related')
    }
    # select largest node
    a = tr1
    b = tr2
    if(Ntip(tr2) > Ntip(tr1)){
        a = tr2
        b = tr1
    }

    d = setRownames(metadata, metadata$id)
    d = d[strip.dx_labs(a$tip.label),]

    # add node to input data
    d$in_node = d$id %in% strip.dx_labs(b$tip.label)
    
    if(type == 'dx'){
         vars = c('(Intercept)','nice_dxIBD')
    }
    if(type == 'full'){
         vars = c('(Intercept)','nice_dxIBD','nice_age','nice_sexF','nice_bmi')
    }
    if(type == 'all'){
         vars = c('IBD','IBD.full','age','sexF','bmi','CD','UC')
    }
    ri = data.frame(var=vars, estiamte= NA, se=NA ,z=NA, pval=NA)
    
    # require minimum n=10 per node
    if((min(table(d$in_node)) > 10) && (min(table(d$nice_dx)) > 10)){
        tryCatch({
            if(type == 'dx'){
                u = summary(glmer(in_node ~ nice_dx + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))$coef
                ri = cbind(rownames(u), u)
            }
            if(type == 'full'){
                u = summary(glmer(in_node ~ nice_dx + nice_age + nice_sex + nice_bmi + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))$coef		    
                ri = cbind(rownames(u), u)
                ri = data.frame(ri)
                ri['nice_dxIBD',1] = 'nice_dxIBD.full'
            }
            if(type == 'all'){
                w = summary(glmer(in_node ~ nice_dx + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))$coef
                ri = cbind(rownames(w), w)
                
                u = summary(glmer(in_node ~ nice_dx + nice_age + nice_sex + nice_bmi + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))$coef		    
                ri = rbind(ri, cbind(rownames(u), u))

                v = summary(glmer(in_node ~ nice_dx2 + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))$coef
                ri = rbind(ri, cbind(rownames(v), v))
                
                ri = data.frame(ri)
                ri['nice_dxIBD.1',1] = 'nice_dxIBD.full'
                
                i = grep('(Intercept)',rownames(ri), value = T, invert = T)
                ri = ri[i, ,drop=F]
                ri[,1] = gsub('nice_|nice_dx|nice_dx2','',ri[,1])
            }
            },
            error = function(...){
                
        }
        )

    }
    colnames(ri) = c('var', 'estimate', 'se', 'z', 'pval')
    ri = setRownames(data.frame(ri), NULL)
    ri[,2:ncol(ri)] = sapply(ri[,2:ncol(ri)], as.numeric)
    return(ri)
}

node.enrichment = function(node, phy, metadata = impute.metadata(), type='all'){
    cbind(node, enrichment.between(node, Ntip(phy)+1, phy, metadata, type))
}
