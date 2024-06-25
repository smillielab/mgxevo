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
parser$add_argument("--tree_dir", help="Path to tree directory")
parser$add_argument("--phylo_enrichment_dir", help="Path to directory containing phylo enrichment results")
parser$add_argument("--output_dir", help="Path to output directory")
args <- parser$parse_args()
tree_dir <- args$tree_dir


simple_select_nodes = function(tree, nodes, min_size=25, other=TRUE){
    
    # get disjoint significant nodes
    nodes.use = c()
    for(ni in nodes){
        flag = 1
	li = extract.clade(tree, ni)$tip.label
	for(nj in nodes.use){
	    lj = extract.clade(tree, nj)$tip.label
	    if(length(intersect(li, lj)) > 0){
	        flag = 0
		break
	    }
	}
	if(flag == 1){
	    nodes.use = c(nodes.use, ni)
	}
    }

    if(other == FALSE){
        return(nodes.use)
    }
    
    # get disjoint other nodes
    n0 = length(tree$tip.label)+1
    nf = n0 + tree$Nnode - 1
    node2size = sort(setNames(sapply(n0:nf, function(a) Ntip(extract.clade(tree, a))), n0:nf), dec=T)
    for(ni in names(node2size)){ni = as.integer(ni)
        if(node2size[[as.character(ni)]] >= min_size){
            flag = 1
	    li = extract.clade(tree, ni)$tip.label
	    for(nj in nodes.use){
	        lj = extract.clade(tree, nj)$tip.label
		if(length(intersect(li, lj)) > 0){
		    flag = 0
		    break
		}
	    }
	    if(flag == 1){
	        nodes.use = c(nodes.use, ni)
	    }
	}
    }
    return(nodes.use)	
}

df = ffread(args$multigene_table)
ph = do.call(rbind, sapply(df$genome, function(gi) fread(sprintf('%s/%s.phylo.enrich.full.txt', args$phylo_enrichment_dir, gi)) ,simplify=F))
ph = ph[test == 'nice_dxIBD' | is.na(test)]
ph[, padj:=p.adjust(pval, 'fdr'), .(name)]
ph = ph[order(padj)][padj < .05][,.(name, node, coef, padj)]
g.use = sort(unique(ph[padj < .05]$name))

# re-run for gene-strain analysis
#ak = ffread(paste0(tree_dir, '/gene_strain.ak.tips.txt'), as.dt=T)
#g.use = sort(unique(ak$name))

tips = sapply(g.use, function(gi){
    tr = read.tree.nice(sprintf('%s/%s.tree',tree_dir, gi))
    nd = ph[padj < .05][name == gi]$node
    #nd = ak[name == gi]$node
    nd = sort(unique(simple_select_nodes(tr, nd, min_size=50, other=FALSE))) # was 25 first run
    sa = sapply(nd, function(a) c(gi, a, paste(extract.clade(tr, a)$tip.label, collapse=';')), simplify=T)
    t(sa)
}, simplify=F)
tips = as.data.table(do.call(rbind, tips))
colnames(tips) = c('name', 'node', 'samples')
tips$node = as.integer(tips$node)

out = merge(ph, tips, by=c('name', 'node'), all=T)[samples != '']
writeLines(g.use, paste0(args$output_dir, '/genomes.use.txt'))
write.table(out, paste0(args$output_dir, '/strainfinder.final_tips.txt'), quote=F, row.names=F, sep='\t')

