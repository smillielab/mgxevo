specific_markers = function(obj, markers, ident.use=NULL, test.use='mast', data.use='log2', tpm=NULL, n.cores=1, ...){

    # Find genes that are DE in ALL pairwise comparisons between (ident.use, obj$ident)
    # Run each test on intersecting genes with keep=TRUE
    # After each test, keep genes ONLY if (pval < .05), (coef > 0), and (log2fc > 0)

    # Check arguments
    markers[, ident := gsub(':.*', '', gsub('ident', '', contrast))]
    markers[, keep := TRUE]
    if(! 'pvalH' %in% colnames(markers) | test.use != 'mast'){stop('test.use != mast')}

    # Get cell identities
    idents = as.character(levels(obj$ident))
    if(is.null(ident.use)){
        ident.use = idents
    }
    remove = setdiff(idents, unique(markers$ident))
    if(length(remove) > 0){print(paste('Removing', paste(remove, sep=', ')))}
    idents = intersect(idents, unique(markers$ident))

    # Pre-calculate TPM
    tpm = get_data(obj, data.use='tpm', tpm=tpm)

    # Pre-calculate data
    data = get_data(obj, data.use=data.use, tpm=tpm)

    # Find marker genes
    for(i in ident.use){

        # Sort comparisons by overlapping genes
        js = setdiff(idents, i)
	js = rev(names(sort(sapply(js, function(j) length(intersect(markers[ident == i, gene], markers[ident == j, gene]))))))

        for(j in js){

	    # get genes to test
	    genes.1 = markers[ident == i & keep == TRUE, gene]
	    genes.2 = markers[ident == j & keep == TRUE, gene]
	    genes.use = intersect(genes.1, genes.2)
	    if(length(genes.use) == 0){next}
	    print(paste(i, j, 'testing', length(genes.use), 'genes'))

	    # get lrt_regex
	    lrt_regex = paste(paste0(unique(markers[ident == i, contrast]), '$'), collapse='|')

	    # run marker test
	    mi = p.find_markers(obj, ident.1=i, ident.2=j, tpm=tpm, data.use=data, genes.use=genes.use, dir='pos', test.use='mast', do.stats=F, min_fc=0, min_alpha=0, lrt_regex=lrt_regex, ...)

	    # update genes to keep
	    if(is.null(mi)){
	        genes.remove = genes.use
	    } else {
	        genes.keep = mi[((pvalD <= .05 & coefD > 0) | (pvalC <= .05 & coefC > 0)) & log2fc > 0, gene]
	        genes.remove = setdiff(genes.use, genes.keep)
	    }

	    print(paste('Removing', paste(genes.remove, collapse=', ')))
	    markers[ident == i & gene %in% genes.remove]$keep = FALSE
	    print(paste('Keeping', paste(markers[ident == i & keep == TRUE, gene], collapse=', ')))
	    print(paste('ident', i, 'found', sum(markers[ident == i, keep]), 'genes'))
	    print(dim(markers))
	}
    }
    markers = markers[ident %in% ident.use,]
    return(markers)
}


reject_markers = function(obj, ident2genes, test.use='mast', data.use='log2', tpm=NULL, n.cores=1, ...){

    # Find genes that are *negatively DE* in at least one pairwise comparison
    # Run each test on genes in ident2genes with reject = FALSE
    # After each test, reject genes ONLY if (pval < .05) & (coef < 0) & (log2fc < 0)

    # Fix input arguments
    res = as.data.table(stack(ident2genes))
    colnames(res) = c('gene', 'ident')
    res[, reject := FALSE]
    res[, reject_ident := 'NA']
    res[, reject_log2fc := 0]
    setkeyv(res, c('ident', 'gene'))

    # Get cell identities
    idents = as.character(levels(obj$ident))
    ident.use = names(ident2genes)

    # Pre-calculate TPM
    tpm = get_data(obj, data.use='tpm', tpm=tpm)

    # Pre-calculate data
    data = get_data(obj, data.use=data.use, tpm=tpm)

    # Find marker genes
    for(i in ident.use){
        for(j in setdiff(idents, i)){

	    # setup marker test
	    genes.use = res[ident == i, gene]
	    lrt_regex = paste0('ident', i, '$')
	    if(length(genes.use) == 0){next}

	    # run marker test
	    mi = p.find_markers(obj, ident.1=i, ident.2=j, tpm=tpm, data.use=data, genes.use=genes.use, dir='neg', test.use='mast', do.stats=F, min_fc=1.01, min_alpha=.005, lrt_regex=lrt_regex, ...)

	    # update genes to keep
	    if(!is.null(mi)){
	        mi = mi[((pvalD <= .05 & coefD < 0) | (pvalC <= .05 & coefC < 0)) & log2fc < 0]
		if(nrow(mi) == 0){next}
		k = data.frame(ident=i, gene=mi[,gene])
		mi = mi[res[k, reject_log2fc] > mi[,log2fc]]
		if(nrow(mi) == 0){next}
		k = data.frame(ident=i, gene=mi[,gene])
		res[k, reject := TRUE]
		res[k, 'reject_ident' := j, with=F]
		res[k, 'reject_log2fc' := mi[,log2fc], with=F]
	    }
	}
    }
    return(res)
}


find_specific_genes = function(markers, anno=NULL, min.n=0, min.same=-Inf, min.next=-Inf, max.pval=.05, keep=FALSE){

    # Identify specific markers as follows:
    # 1) Do minimal filtering of the marker list to calculate "next stats"
    # 2) Filter by "next stats"
    # 3) Finally, filter by p-value to make sure each marker is significant
    # Optionally, use "keep" column to filter while calculating "next stats"

    # Process input data
    markers[, specific := FALSE]
    markers = markers[n >= min.n]

    # Calculate next stats
    f = ifelse(keep == FALSE, identity, function(m){m[keep == TRUE]})
    markers = next_stats(markers, anno=anno, filter=f)

    # Filter by p-value
    i = quote((coefD >= 0 & pvalD <= max.pval) | (coefC >= 0 & pvalC <= max.pval))
    j = quote((same_log2fc >= min.same) & (next_log2fc >= min.next))
    markers[eval(i) & eval(j), specific := TRUE]
    if(keep == TRUE & !is.infinite(min.next)){markers[keep == FALSE, specific := FALSE]}

    markers

}


