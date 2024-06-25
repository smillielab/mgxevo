require(data.table)
source('/broad/smillie-data/code/csmillie/single_cell/tpm.r')


fast_gsea = function(ranks, pathways=NULL, bg=NULL, do.bg=FALSE, regex='BP|MF|CC|hallmark|canonical', minSize=5, maxSize=500, nperm=1000, gseaParam=1){

    require(fgsea)
    hpath = '/broad/smillie-data/code/csmillie/db/human.gsea_pathways.rds'
    mpath = '/broad/smillie-data/code/csmillie/db/mouse.gsea_pathways.rds'

    # inputs:
    # ranks = list(gene1 = pval, gene2 = pval)
    # pathways = list(db1 = list(sig1 = genes, sig2 = genes, ...), db2 = list(sig1 = genes, sig2 = genes, ...))
    # bg = background set of genes (remove these from ranks and pathways)

    # notes:
    # this function should be run on ranks from full gene universe, e.g.
    # markers = markers[pvalH < .05]
    # ranks = structure(markers$log2fc, names=markers$gene)

    # convert markers to ranks
    if('gene' %in% colnames(ranks) & 'coefD' %in% colnames(ranks)){
        ranks = ranks[pvalD < .05]
	ranks = structure(ranks$coefD, names=ranks$gene)
    }

    # load default pathways
    if(is.null(pathways)){
        org = predict_organism(names(ranks)[1:min(length(ranks), 10)])
	if(org == 'human'){pathways = readRDS(hpath)} else {pathways = readRDS(mpath)}
    } else if(!is.list(pathways)){
        if(length(pathways) == 1){
	    if(pathways == 'kegg'){
	        pathways = list(KEGG=load_kegg(do.flatten=T))
		regex = ''
	    }
	} else {
	    pathways = list(Pathway=pathways)
	}
    }
    pathways = pathways[grep(regex, names(pathways))]
    print(names(pathways))

    # filter background set
    if(is.null(bg)){bg = names(ranks)}
    bg = intersect(bg, names(ranks))
    if(do.bg == TRUE){
        ranks = ranks[bg]
	pathways = lapply(pathways, function(a) lapply(a, function(b) intersect(bg, b)))
    }

    gsea = list()

    # run gsea on each set of pathways
    for(name in names(pathways)){

        # filter pathways by size
	pathway = pathways[[name]]
	sizes = sapply(pathway, length)
	pathway = pathway[minSize <= sizes & sizes <= maxSize]
	print(paste(name, 'testing', length(pathway), 'gene sets'))

	# run gsea and sort by NES
        res = fgsea(pathways=pathway, stats=ranks, nperm=nperm, minSize=minSize, maxSize=maxSize, gseaParam=gseaParam)
	res = res[order(-1*res$NES),]
	gsea[[name]] = res
    }

    return(gsea)
}


fast_enrich = function(genes, regex='GO_.*2017$|KEGG.*2016|Reactome.*2016|Panther_2016', collapse=FALSE){

    require(enrichR)
    require(data.table)

    # For full list of databases: listEnrichrDbs()
    # Computes enrichment with Fisher exact test
    # Also uses random gene sets to correct the Fisher test

    # Select databases to use
    dbs = grep(regex, listEnrichrDbs()[,3], value=T)

    # Run enrichment test
    res = enrichr(genes, dbs)

    # Fix each list
    res = sapply(res, function(a) {

        # Sort by adjusted p-value
        a = as.data.table(a)[,.(Term, Overlap, P.value, Adjusted.P.value, Genes)][order(Adjusted.P.value)]

	# Get overlap statistics
	a[, NM := as.integer(gsub('/.*', '', Overlap))]
	a[, NQ := length(genes)]
	a[, NT := as.integer(gsub('.*/', '', Overlap))]

    }, simplify=F)

    # Return results
    if(collapse == TRUE){res = do.call(rbind, res)[order(P.value)]}
    res
}


gsea_heatmap = function(gsea, markers=NULL, colors=NULL, max_pval=.05, max_terms=50, max_genes=200, min_genes_per_term=2, max_genes_per_term=20, font_size=8, fix_names=TRUE, max_overlap=.5,
                        xstag=FALSE, ystag=FALSE, xsec=FALSE, ysec=FALSE, auto_condense='cols', out=NULL, nrow='auto', ncol='auto', dir='pos', ...){
    require(tibble)

    # Plot heatmap of GSEA results (terms x genes)
    # automatically select terms and genes within specified bounds
    # optionally color using named list of genes

    # Get data to plot
    if('Term' %in% colnames(gsea) & 'Genes' %in% colnames(gsea)){
        terms = gsea$Term
	genes = gsea$Genes
    }

    if('pathway' %in% colnames(gsea) & 'leadingEdge' %in% colnames(gsea)){
        gsea = gsea[pval <= max_pval,]
	if(nrow(gsea) == 0){return(NULL)}
	if(dir == 'pos'){gsea = gsea[NES >= 0]}
	if(dir == 'neg'){gsea = gsea[NES <= 0]}
	terms = gsea$pathway
	genes = gsea$leadingEdge
	pvals = -log10(gsea$pval)
    }

    # Split genes
    if(! is.list(genes)){genes = strsplit(genes, ',|;| ')}

    # Remove redundant terms
    remove = sapply(2:length(genes), function(i){
        any(sapply(1:(i-1), function(j) {
	    a = length(intersect(genes[[i]], genes[[j]]))
	    b = length(genes[[i]])
	    a/b > max_overlap
	}))
    })
    terms = terms[c(TRUE, !remove)]
    genes = genes[c(TRUE, !remove)]
    pvals = pvals[c(TRUE, !remove)]

    # Filter by minimum genes per term
    i = sapply(genes, length) >= min_genes_per_term
    terms = terms[i]
    genes = genes[i]
    pvals = pvals[i]
    if(length(terms) <= 1 | length(genes) <= 1){return(NULL)}

    # Select max_genes_per_term from each gene list
    genes = sapply(genes, function(a) a[1:min(length(a), max_genes_per_term)])

    # Filter by maximum total genes
    for(i in 1:length(genes)){if(length(unique(unlist(genes[1:i]))) >= max_genes){break}}
    terms = terms[1:(i-1)]
    genes = genes[1:(i-1)]
    pvals = pvals[1:(i-1)]
    if(length(terms) <= 1 | length(genes) <= 1){return(NULL)}

    # Filter by maximum total terms
    if(length(terms) >= max_terms){
        print(paste('Dropping', max_terms - length(terms), 'terms'))
        terms = terms[1:max_terms]
	genes = genes[1:max_terms]
	pvals = pvals[1:max_terms]
    }

    # Get all genes
    all_genes = sort(unique(unlist(genes)))

    # Fix names
    if(fix_names == TRUE){
        terms = fix_gsea_names(terms)
    }

    # Get colors
    if(is.null(colors) & !is.null(markers)){colors = structure(markers$coefD, names=markers$gene); legend.title='coefD'}
    if(is.null(colors)){colors = structure(rep(1, length(all_genes)), names=all_genes); legend.title='Presence'}

    # Incidence matrix
    x = t(sapply(genes, function(gi) all_genes %in% gi) * colors[all_genes])
    rownames(x) = terms
    colnames(x) = all_genes
    if(nrow(x) <= 1 | ncol(x) <= 1){return(NULL)}

    # Heatmap [terms x genes]
    if(ncol == 'auto' & nrow == 'auto'){
        ri = nrow(x)/max(nrow(x), ncol(x))
	ci = ncol(x)/max(nrow(x), ncol(x))
	nrow = max(2.5*ri, 1.25)
	ncol = max(2.5*ci, 1.25)
	print(c(nrow, ncol))
    }

    # Auto condense labels
    if(auto_condense %in% c('rows', 'both')){
        ysec = FALSE
	if(nrow(x) >= max_terms/2){ystag = TRUE}
    }
    if(auto_condense %in% c('cols', 'both')){
        xsec = FALSE
	if(ncol(x) >= max_genes/2){xstag = TRUE}
    }

    p1 = ggheatmap(x, pal='Blues', font_size=font_size, Rowv='hclust', Colv='hclust', xstag=xstag, ystag=ystag, xsec=xsec, ysec=ysec, ...)

    # P-value barplot
    x = data.frame(Term=factor(terms, levels=levels(p1$data$row)), P.value=pvals)
    p2 = ggplot(x, aes(x=Term, y=P.value)) +
         geom_bar(stat='identity', width=.6, fill='#6BAED6') +
	 coord_flip() + xlab('') + ylab('-log10(P-value)') + theme_cowplot(font_size=font_size) +
         theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), panel.grid.major.x=element_line(colour='#000000', size=.25, linetype='dotted')) +
	 scale_y_continuous(position='top', breaks=function(a) seq(floor(min(a)), floor(max(a))), labels=function(a) floor(a))

    # Merge plots
    p3 = get_legend(p1)
    p1 = p1 + theme(legend.position='none')
    p = plot_grid(plot_grid(p1, p2, nrow=1, rel_widths=c(.925, .075), align='h'), p3, nrow=1, rel_widths=c(.95, .05))

    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }
    p
}


gsea_heatmaps = function(gsea, markers=NULL, colors=NULL, prefix, ...){
    outdir = dirname(prefix)
    if(!dir.exists(outdir)){dir.create(outdir, recursive=TRUE)}
    for(name in names(gsea)){
        out = gsub('\\/\\.', '/', paste(prefix, name, 'gsea.png', sep='.'))
	print(paste('Writing', name, 'heatmap to', out))
        gsea_heatmap(gsea=gsea[[name]], markers=markers, colors=colors, out=out, ...)
    }
}


fix_gsea_names = function(terms, nlen=45){
    terms = gsub('GO[^\\)]*;', '', terms)
    terms = gsub('positive', 'pos', terms)
    terms = gsub('negative', 'neg', terms)
    terms = gsub('regulation', 'reg', terms)
    terms = gsub('response', 'resp', terms)
    terms = gsub('signaling', 'sig', terms)
    terms = gsub('interaction', 'ix', terms)
    terms = gsub('DOWN', 'down', gsub('UP', 'up', terms))
    terms = make.unique(substr(terms, 1, nlen))
    terms
}


gsealist_heatmap = function(glist, regex='', name='BP', n=3, n.plot=Inf, max_pval=1, max_padj=1, plot.max_pval=1, plot.max_padj=1, fix_names=TRUE, nlen=50, replace_na=0,
                            rcut=NULL, ret.data=F, ...){

    # Filter by p-value and sort by NES
    glist = glist[grep(regex, names(glist))]
    glist = lapply(glist, function(a) a[[name]][order(-1*NES)])
    glist = glist[lapply(glist, nrow) > 0]

    # Select [n] terms from each gsea
    terms = sort(unique(unlist(as.vector(sapply(glist, function(gsea) gsea[pval <= max_pval & padj <= max_padj][NES > 0][1:min(n, nrow(gsea))]$pathway)))))

    # Get [terms x NES] matrix
    data = sapply(glist, function(gsea){
        gsea = gsea[pval <= plot.max_pval & padj <= plot.max_padj][1:min(n.plot, nrow(gsea))]
	setkey(gsea, pathway)
	gsea[terms, NES]
    })

    # Get [terms x padj] matrix
    pvals = sapply(glist, function(gsea){
        gsea = gsea[pval <= plot.max_pval & padj <= plot.max_padj][1:min(n.plot, nrow(gsea))]
	setkey(gsea, pathway)
	gsea[terms, pval]
    })

    # Fix rownames and NAs
    if(fix_names == TRUE){terms = fix_gsea_names(terms, nlen=nlen)}
    rownames(data) = rownames(pvals) = terms
    data[is.na(data)] = replace_na
    pvals[is.na(pvals)] = 1

    # Remove redundant terms with rcut
    if(!is.null(rcut)){
        keep = rep(TRUE, nrow(data))
    	dist = cor(t(data))
    	for(i in 1:nrow(data)){
            if(keep[[i]] == FALSE){next}
            j = which(dist[i,] > rcut)
            if(length(j) <= 1){next}
            k = apply(data[j,], 1, quantile, .8)
            j = j[k != max(k)]
            keep[j] = FALSE
    	}
    	data = data[keep,]
    }

    if(ret.data == TRUE){list(data=data, pvals=pvals)} else {
        ggheatmap(data, pvals=pvals, max_pval=.05, ...)
    }
}


score_signatures2 = function(obj, signatures, group_by=NULL, cells.use=NULL){

    # Fix input arguments
    data = obj$data
    if(is.null(cells.use)){cells.use = colnames(data)}
    if(!is.null(group_by)){
        group_by = as.factor(group_by[cells.use])
	totals = table(group_by)
    }

    # Calculate signature scores
    sapply(signatures, function(genes.use){
        genes.use = intersect(genes.use, rownames(data))
	x = as.matrix(t(data[genes.use, cells.use]))
	if(!is.null(group_by)){
	    x = rowsum(x, group_by)
	    if(any(rownames(x) != names(totals))){stop()}
	    x = x/as.vector(totals)
	}
	x = x[,apply(x, 2, sd) != 0]
	rowMeans(x)
   })
}

diff_signatures = function(obj, cells.1, cells.2, signatures, group_by=NULL){

    # Fix input arguments
    data = obj$data
    if(is.null(group_by)){group_by = obj$ident}
    if(is.null(names(group_by))){stop()}
    groups.1 = as.factor(group_by[as.character(cells.1)])
    groups.2 = as.factor(group_by[as.character(cells.2)])
    totals.1 = table(groups.1)
    totals.2 = table(groups.2)

    sapply(signatures, function(genes.use){
        genes.use = intersect(genes.use, rownames(data))
	x1 = rowsum(as.matrix(t(data[genes.use, cells.1])), groups.1)
	x2 = rowsum(as.matrix(t(data[genes.use, cells.2])), groups.2)
	if(any(rownames(x1) != names(totals.1)) | any(rownames(x2) != names(totals.2))){stop()}
	x1 = x1/as.vector(totals.1)
	x2 = x2/as.vector(totals.2)
	diff = x1 - x2
	rowMeans(diff)
    })
}


pair_enrich = function(pairs, gene_sets, universe='/broad/smillie-data/code/csmillie/db/hg19_genes.txt', require_ix=TRUE){
    require(data.table)

    # Calculate enrichment scores for all gene pairs (e.g. ligand-receptor interactions) and gene sets (e.g. GO terms)
    # This function uses the minimum overlap found in each column of the gene pair
    # If require_ix, then both genes in the gene pair must be found in the gene set

    # Read input data
    universe = readLines(universe)
    pairs = as.data.table(pairs)
    colnames(pairs) = c('a', 'b')

    # Select gene universe
    pairs = pairs[a %in% universe & b %in% universe]
    gene_sets = sapply(gene_sets, function(a) intersect(a, universe), simplify=F)

    # Pre-calculate sizes
    nu = length(unique(universe))
    np = length(unique(unlist(pairs)))

    # Fisher test on gene pairs
    res = sapply(names(gene_sets), function(name){
        gene_set = gene_sets[[name]]

	if(require_ix == TRUE){
	    p = pairs[a %in% gene_set & b %in% gene_set]
	} else {
	    p = pairs
	}

	# Calculate intersections
	a = unique(intersect(p$a, gene_set))
	b = unique(intersect(p$b, gene_set))

	# (2x2) contingency table
	u = 2*min(length(a), length(b))
	v = length(gene_set) - u
	w = np - u
	x = nu - u - v - w

	# Fisher test
	if(u < 2){return(NULL)}
	data.frame(Term=name, P.value=fisher.test(matrix(c(u,v,w,x), nrow=2, ncol=2))$p.value, Ligand=paste(a, collapse=','), Receptor=paste(b, collapse=','), N=u)
    })

    as.data.table(do.call(rbind, res))[order(P.value)]
}


gsea.fisher = function(gene_set1, gene_set2, universe='/broad/smillie-data/code/csmillie/db/hg19_genes.txt', n.cores=1){

    # GSEA with fisher test
    # gene_set1 = list of genes
    # gene_set2 = list of genes
    # return matrix of p-values

    # convert gene sets to list
    if(typeof(gene_set1) != 'list'){gene_set1 = list(gene_set1)}
    if(typeof(gene_set2) != 'list'){gene_set2 = list(gene_set2)}

    # fix names
    if(is.null(names(gene_set1))){names(gene_set1) = paste0('query_', 1:length(gene_set1))}
    if(is.null(names(gene_set2))){names(gene_set2) = paste0('target_', 1:length(gene_set2))}

    # pairwise fisher tests
    if(length(universe) == 1){all_genes = readLines(universe)} else {all_genes = universe}
    m = run_parallel(
        foreach(i=names(gene_set1), .combine=rbind) %:% foreach(j=names(gene_set2), .combine=rbind) %dopar% {
	    gi = gene_set1[[i]]
	    gj = gene_set2[[j]]
	    u = factor(all_genes %in% unlist(gi), levels=c(FALSE, TRUE))
	    v = factor(all_genes %in% unlist(gj), levels=c(FALSE, TRUE))
	    q = fisher.test(table(u,v))
	    o = intersect(gi, gj)
	    d = data.frame(query=i, target=j, pval=q$p.value, odds_ratio=q$estimate, overlap=paste(o, collapse=','))
	    d
	},
	n.cores=n.cores
    )
    m = as.data.table(m)
    m[, padj := p.adjust(pval, 'fdr')]
    m = m[order(pval)]
    return(m)
}


gsea.mast = function(data, covariates, formula=NULL, lrt_regex=TRUE, gsea.boot=100, n.cores=1){

    load_mast()
    options(mc.cores=n.cores)

    # Make single cell assay object
    fdata = data.frame(matrix(rep(1, nrow(data))))
    covariates = as.data.frame(covariates)
    sca = MAST::FromMatrix(as.matrix(data), covariates, fdata)

    # Fit MAST hurdle model
    if(is.null(formula)){
        formula = paste('~' , paste(colnames(covariates), collapse=' + '))
    }
    zlm.obj = zlm(as.formula(formula), sca)

    # Calculate bootstraps
    boot = bootVcov1(zlm.obj, gsea.boot)

    # Get GO gene list
    genes = read.table('/broad/smillie-data/code/csmillie/db/gopca/go_annotations_human.tsv', sep='\t', stringsAsFactors=F, row.names=1)
    sets = structure(strsplit(genes[,4], ','), names=rownames(genes))
    sets = lapply(sets, function(a){b = as.integer(na.omit(match(a, rownames(zlm.obj@coefC))))})
    sets = sets[sapply(sets, length) >= 5]

    # Get hypothesis columns
    names = colnames(zlm.obj@coefC)[grep(paste(lrt_regex, collapse='|'), colnames(zlm.obj@coefC))]

    # Perform GSEA
    res = lapply(names, function(a){
        gsea = summary(gseaAfterBoot(zlm.obj, boot, sets, CoefficientHypothesis(a)))
	gsea$name = genes[as.character(gsea$set), 3]
	gsea$genes = genes[as.character(gsea$set), 4]
	gsea$contrast = a
	return(gsea)
    })
    res = do.call(rbind, res)

    options(mc.cores=1)
    return(res)
}


go_genes = function(genes, universe, ontology='BP'){

    require(topGO)
    GO2genes = readMappings(file='~/aviv/db/gopca/go.test.txt', sep='\t')
    gene_list = as.numeric(universe %in% genes)
    names(gene_list) = universe
    gene_list = factor(gene_list)

    # Run topGO tests
    GOdata = new('topGOdata', ontology=ontology, allGenes=gene_list, annot=annFUN.GO2genes, GO2genes=GO2genes)
    res = runTest(GOdata, algorithm='classic', statistic='fisher')
    res = GenTable(GOdata, res, topNodes=min(1000, length(res@score)))
    res = res[res$result1 <= .05,]
    return(res)
}


go_markers = function(m, top=NULL, pval=NULL, auc=NULL, ontology='BP', n.cores=1){

    require(topGO)
    require(naturalsort)
    GO2genes = readMappings(file='~/aviv/db/gopca/go.test.txt', sep='\t')

    # Ontology can be: BP (biological process), MF (molecular function), or CC (cellular component)

    # Get column names
    if('contrast' %in% colnames(m)){
	score = 'pval'
	cluster = 'contrast'
	m = m[m$log2fc > 0,]
    } else {
	score = 'auc'
	cluster = 'cluster'
	m = m[m$auc > .5,]
    }

    # Test each cluster
    all_genes = unique(m$gene)
    clusters = naturalsort(as.character(unique(m[,cluster])))

    go_terms = run_parallel(foreach(a=clusters, .packages=c('topGO')) %dopar% {

        # Filter marker genes
        mi = m[m[,cluster] == a,]
	if(!is.null(top)){
	    print('Using top genes')
	    mi = mi[1:top,]
	}
	if(!is.null(pval)){
	    print('Filtering by pval')
	    mi = mi[mi$pval <= pval,]
	}
	if(!is.null(auc)){
	    print('Filtering by AUC')
	    mi = mi[mi$auc <= auc,]
	}

	# Construct gene list
	gene_list = as.numeric(all_genes %in% mi$gene)
	names(gene_list) = all_genes
	gene_list = factor(gene_list)

	# Run topGO tests
	GOdata = new('topGOdata', ontology=ontology, allGenes=gene_list, annot=annFUN.GO2genes, GO2genes=GO2genes)
	res = runTest(GOdata, algorithm='classic', statistic='ks')
	res = GenTable(GOdata, res, topNodes=min(1000, length(res@score)))
	res = res[res$result1 <= .05,]
	return(res)
    }, n.cores=n.cores)

    names(go_terms) = clusters
    return(go_terms)
}

load_kegg = function(names.use=NULL, names.rmv=NULL, do.flatten=FALSE, no_spaces=FALSE, regex=NULL){
    kegg = readRDS('/broad/smillie-data/db/kegg/3.kegg.human.rds')

    kegg = sapply(names(kegg), function(A)
               sapply(names(kegg[[A]]), function(B)
	           sapply(names(kegg[[A]][[B]]), function(C) {
		       if(!is.null(regex)){if((! grepl(regex, A) & ! grepl(regex, B)) & ! grepl(regex, C)){return(NULL)}}
		       if(!is.null(names.use)){if(!(A %in% names.use | B %in% names.use | C %in% names.use)){return(NULL)}}
		       if(!is.null(names.rmv)){if(A %in% names.rmv | B %in% names.rmv | C %in% names.rmv){return(NULL)}}
		       as.character(na.omit(kegg[[A]][[B]][[C]]))
		   }, simplify=F),
	       simplify=F),
	   simplify=F)

    if(no_spaces == TRUE){
        for(A in names(kegg)){
	    for(B in names(kegg[[A]])){
	        names(kegg[[A]][[B]]) = gsub(' ', '_', names(kegg[[A]][[B]]))
	    }
	    names(kegg[[A]]) = gsub(' ', '_', names(kegg[[A]]))
	}
	names(kegg) = gsub(' ', '_', names(kegg))
    }

    if(do.flatten == TRUE){
        kegg = unlist(unlist(kegg, recursive=FALSE), recursive=FALSE)
	kegg = kegg[!sapply(kegg, is.null)]
    }
    kegg
}
