
expand.grid.df = function(...) Reduce(function(...) merge(..., by=NULL), list(...))

smart_shuffle = function(data, cols.shuf, cols.test=NULL, max_tries=10, max_iter=100){

    # Shuffle columns while avoiding duplicate pairs
    # data = [m x n] matrix
    # cols.shuf = list of columns to shuffle
    # cols.test = list of columns to test

    # Fix input data
    data = as.data.frame(data)
    if(is.null(cols.test)){cols.test = colnames(data)}
    if(! all(cols.shuf %in% cols.test)){stop('error: cols.shuf not in cols.test')}

    # For each attempt, shuffle columns
    for(a in 1:max_tries){
	data[, cols.shuf] = data[sample(1:nrow(data)), cols.shuf]

	# Iteratively fix duplicates
	for(b in 1:max_iter){
	    new_data = data[, cols.shuf, drop=F]
	    i = duplicated(data[,cols.test])
	    if(sum(i) == 0){return(as.data.table(data))}
	    j = sample(which(!i), sum(i))
	    new_data[i,] = data[j, cols.shuf]
	    new_data[j,] = data[i, cols.shuf]
	    data[, cols.shuf] = new_data
	}
    }
    as.data.table(data)
}

make_edges = function(data, diff=NULL, weights=NULL, ix=NULL, do.intersect=FALSE, full_complex=FALSE){
    require(data.table)

    # Build cell-cell interaction network from markers and gene pairs
    #
    # Input arguments:
    # - data = matrix(1=ident, 2=gene) = cell type markers
    # - diff = matrix(1=ident, 2=gene) = DE genes (optional)
    # - weights = matrix(1=ident, 2=gene, 3=weight)
    # - ix = matrix(1=gene, 2=gene) = interaction matrix
    # - permute = permute gene list?
    # - do.intersect = intersect diff with data, i.e.
    #   require DE genes to also be cell type markers
    # - returns edgelist for data or diff

    # Fix inputs
    data = as.data.table(data[,c('ident', 'gene')])
    data[, ct := 1]

    # Merge markers with DE genes
    if(!is.null(diff)){
        diff = as.data.table(diff[,c('ident', 'gene')])
	diff[, de := 1]
	data = merge(data, diff, by=c('ident', 'gene'), all=TRUE)
    } else {
        data[, de := 0]
    }
    data = setkeyv(data, c('ident', 'gene'))

    # Remove DE genes that are not cell type markers
    if(do.intersect == TRUE){data = data[ct == 1]}

    # Read interactions
    if(is.null(ix)){
        ix = read.table('~/code/single_cell/PairsLigRec.txt', sep='\t', header=T, stringsAsFactors=F, comment.char='', quote='')
	ix = ix[ix[,'Pair.Evidence'] == 'literature supported',]
	ix = ix[,c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol')]
	colnames(ix) = c('lig', 'rec')
    }
    ix = as.data.table(sapply(ix, as.character))
    colnames(ix) = c('lig', 'rec')
    ix$lig = strsplit(ix$lig, ',')
    ix$rec = strsplit(ix$rec, ',')

    if(full_complex == FALSE){
        ix = as.data.table(unique(do.call(rbind, apply(ix, 1, function(a) expand.grid(a[[1]], a[[2]])))))
    	ix = as.data.table(sapply(ix, as.character))
        colnames(ix) = c('lig', 'rec')
    }

    if(is.null(diff)){
        i = apply(ix, 1, function(a) any(a[[1]] %in% data$gene) & any(a[[2]] %in% data$gene))
        ix = ix[i]
    } else {
        i = apply(ix, 1, function(a) any(a[[1]] %in% data$gene) & any(a[[2]] %in% diff$gene))
	j = apply(ix, 1, function(a) any(a[[1]] %in% diff$gene) & any(a[[2]] %in% data$gene))
	ix = ix[i | j]
    }
    genes.use = sort(unique(c(unlist(ix$lig), unlist(ix$rec))))
    data = data[gene %in% genes.use]

    # Initialize network
    nodes = sort(unique(data$ident))
    edges = c()
    graph = matrix(0, nrow=length(nodes), ncol=length(nodes))
    rownames(graph) = colnames(graph) = nodes
    n = list(edges=edges, graph=graph)
    d = list(edges=edges, graph=graph)

    # Iterate over ligand-receptor pairs
    for(i in 1:nrow(ix)){

        # Ligands and receptors
	l = ix[i, lig]
	r = ix[i, rec]

	# Baseline edges
	n.edges = c()
	n.l = data[ct == TRUE][, all(l %in% gene), ident][V1 == TRUE, ident]
	n.r = data[ct == TRUE][, all(r %in% gene), ident][V1 == TRUE, ident]
	if(length(n.l) > 0 & length(n.r) > 0){
	    n.edges = as.matrix(expand.grid(n.l, n.r))
	}

	# DE edges
	d.edges = c()
	d.l = data[de == TRUE][, all(l %in% gene), ident][V1 == TRUE, ident]
	D.l = unique(c(n.l, d.l))
	d.r = data[de == TRUE][, all(r %in% gene), ident][V1 == TRUE, ident]
	D.r = unique(c(n.r, d.r))

	if(length(d.l) > 0 & length(D.r) > 0){
	    d.edges = rbind(d.edges, as.matrix(expand.grid(d.l, D.r)))
	}
	if(length(D.l) > 0 & length(d.l) > 0){
	    d.edges = rbind(d.edges, as.matrix(expand.grid(D.l, d.r)))
	}
	d.edges = unique(d.edges)

	# Fix names
	if(!is.null(nrow(n.edges))){if(nrow(n.edges) > 0){
	    n.edges = cbind(n.edges, paste(l, collapse=','), paste(r, collapse=','))
	    colnames(n.edges) = c('lcell', 'rcell', 'lig', 'rec')
	    n$edges = rbind(n$edges, n.edges)
	}}
	if(!is.null(nrow(d.edges))){if(nrow(d.edges) > 0){
	    d.edges = cbind(d.edges, paste(l, collapse=','), paste(r, collapse=','))
	    colnames(d.edges) = c('lcell', 'rcell', 'lig', 'rec')
	    d$edges = rbind(d$edges, d.edges)
	}}
    }

    # Build graph from edgelist
    if(is.null(diff)){
        edges = as.data.table(n$edges)
    } else {
        edges = as.data.table(d$edges)
    }

    return(edges)
}


permute_edges = function(edges, perm.col='ident', groups=NULL){
    require(data.table)

    # Get column to permute
    if(perm.col == 'ident'){
        perm.l = 'lcell'
	perm.r = 'rcell'
    } else {
        perm.l = 'lig'
	perm.r = 'rec'
    }

    # Do shuffle within groups
    if(is.null(groups)){groups = list(unique(c(edges[[perm.l]], edges[[perm.r]])))}

    # Permute ligands
    map = unique(edges[,.(lcell, lig)])
    new = map
    for(g in groups){
        i = map$lcell %in% g
	new[i] = smart_shuffle(map[i], perm.l, c('lcell', 'lig'))
    }
    #new = smart_shuffle(map, perm.l, c('lcell', 'lig'))
    map[,new_lcell := new$lcell]
    setkeyv(map, c('lcell', 'lig'))
    edges$lcell = map[edges[,.(lcell, lig)], new_lcell]

    # Permute receptors
    map = unique(edges[,.(rcell, rec)])
    new = map
    for(g in groups){
        i = map$rcell %in% g
	new[i] = smart_shuffle(map[i], perm.r, c('rcell', 'rec'))
    }
    #new = smart_shuffle(map, perm.r, c('rcell', 'rec'))
    map[,new_rcell := new$rcell]
    setkeyv(map, c('rcell', 'rec'))
    edges$rcell = map[edges[,.(rcell, rec)], new_rcell]

    # Return edges
    return(edges)
}


edges2graph = function(edges, filter_fxn=identity, symm=TRUE, unique=TRUE, node_order=NULL, permute=FALSE, perm.col='ident', groups=NULL){
    require(data.table)

    # Convert edgelist to adjacency matrix
    # edges = list(1=lcell, 2=rcell, 3=lig, 4=rec, 5=lweight, 6=rweight, 7+=columns for edge filtering)
    # filter = function for filtering edges
    # symm = make graph symmetric by: x = x + t(x)
    # unique = count unique ligands and receptors
    # groups = list(c(ident1, ident2, ident3), c(ident4, ident5), ...) shuffles ligands within cell groups

    # Check formatting
    edges = as.data.table(edges)
    if(! all(c('lcell', 'rcell', 'lig', 'rec', 'lweight', 'rweight') %in% colnames(edges))){
        stop('colnames(edges) != lcell, rcell, lig, rec, lweight, rweight')
    }

    # Permute ligands and receptors
    if(permute == TRUE){
        edges = permute_edges(edges, perm.col=perm.col, groups=groups)
    }

    # Filter edges
    edges = filter_fxn(edges)

    # Get nodes and genes
    if(is.null(node_order)){
        nodes = sort(unique(c(edges$lcell, edges$rcell)))
    } else {
        nodes = node_order
    }
    genes = sort(unique(c(edges$lig, edges$rec)))

    # Initialize adjacency matrix
    graph = as.data.frame(matrix(0, nrow=length(nodes), ncol=length(nodes)))
    rownames(graph) = colnames(graph) = nodes

    if(unique == TRUE){

        # Get unique ligands and receptors
    	lig = unique(edges[,.(lcell, rcell, lig, lweight)])
    	rec = unique(edges[,.(lcell, rcell, rec, rweight)])

    	# Aggregate edge weights
    	lig = lig[, .(weight=sum(lweight)), .(lcell, rcell)]
    	rec = rec[, .(weight=sum(rweight)), .(lcell, rcell)]

    	# Convert to matrix
        lig = data.frame(spread(lig, rcell, weight), row.names=1)
        rec = data.frame(spread(rec, rcell, weight), row.names=1)
        lig[is.na(lig)] = 0
        rec[is.na(rec)] = 0

        # Align data and combine
        x = y = graph
        x[rownames(lig), colnames(lig)] = lig
        y[rownames(rec), colnames(rec)] = rec
        graph = as.matrix(x + y)

    } else {

        # Aggregate edge weights
	edges = edges[, .(weight=sum(lweight + rweight)), .(lcell, rcell)]

	# Convert to matrix
	edges = data.frame(spread(edges, rcell, weight), row.names=1)
	edges[is.na(edges)] = 0

	# Align data and combine
	graph[rownames(edges), colnames(edges)] = edges
    }

    # Make symmetric
    if(symm == TRUE){
        graph = graph + t(graph)
    }

    return(graph)
}

matrix_pvals = function(graph, shuf, type='gt'){

    # Calculate empirical p-values
    # true = [m x n] matrix
    # shuf = list([m x n] matrices)
    # type 'gt' = test true data > shuffled data
    # returns [m x n] matrix of p-values

    if(type == 'lt'){
        Reduce('+', sapply(shuf, function(a) graph >= a, simplify=F))/length(shuf)
    }
    if(type == 'gt'){
        Reduce('+', sapply(shuf, function(a) graph <= a, simplify=F))/length(shuf)
    }
}

matrix_mean = function(graphs){Reduce('+', graphs)/length(graphs)}

matrix_qvals = function(graph, shuf, type='gt', n=100, pvals=c(0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2)){

    true_pvals = matrix_pvals(graph, shuf)
    true_pvals = true_pvals[lower.tri(true_pvals, diag=T)]

    shuf_pvals = as.vector(sapply(1:min(n, length(shuf)), function(i){print(i)
        p = matrix_pvals(shuf[[i]], shuf[setdiff(1:length(shuf), i)])
        p = p[lower.tri(p, diag=T)]
        as.vector(p)
    }))
    num = sapply(pvals, function(p) sum(true_pvals <= p))
    fdr = sapply(pvals, function(p) mean(shuf_pvals <= p)/mean(true_pvals <= p))
    return(list(pvals=pvals, num=num, fdr=fdr))
}

fast_ix = function(data){

    # data = matrix(1=ident, 2=gene) = cell type markers
    edges = make_edges(data)
    graph = edges2graph(edges)
    shuffled_graphs = lapply(1:100, function(i) edges2graph(edges, permute=T))
    pvals = matrix_pvals(graph, shuffled_graphs)
    d = pvals < .05
    nice_network(d)
}
