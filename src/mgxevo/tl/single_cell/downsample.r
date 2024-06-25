require(rsvd)
require(cccd)
require(expm)


num_cells_per_group = function(groups, total_cells=NULL, cells_per_group=NULL){

    num_cells = sort(table(groups))

    if(!is.null(cells_per_group)){
        num_cells[num_cells > cells_per_group] = cells_per_group
    } else {
        n = sort(table(groups))
	if(length(n) == 1){
	    num_cells = total_cells
	    names(num_cells) = names(n)
	} else {
	    u = c(0, cumsum(n)[1:(length(n)-1)])
	    i = (total_cells - u)/seq(length(n), 1, -1) < n
	    if(sum(i) > 0){
	        num_cells[i] = as.integer(ceiling((total_cells - sum(n[!i]))/sum(i)))
	    }
	}
    }

    num_cells
}

resample = function(x,...){if(length(x)==1) x else sample(x,...)}

simple_downsample = function(cells, groups, ngene=NULL, total_cells=NULL, cells_per_group=NULL, replace=FALSE){

    # Downsample cells evenly across groups
    # Option: provide ngene to select HQ cells

    # Set ngene (default = 1)
    if(is.null(ngene)){
        ngene = structure(rep(1, length(cells)), names=cells)
    }
    if(is.null(names(ngene))){names(ngene) = cells}

    # Calculate group sizes
    groups = as.factor(groups)
    num_cells_per_group = num_cells_per_group(groups=groups, total_cells=total_cells, cells_per_group=cells_per_group)

    # Downsample cells within each group
    ds.cells = sapply(levels(groups), function(a){

        # Shuffle cells within group
        cells = resample(cells[groups == a], replace=replace)

	# Select by highest ngene
	cells[order(ngene[cells], decreasing=T)[1:num_cells_per_group[[a]]]]
    })
    ds.cells = as.character(na.omit(unname(unlist(ds.cells))))

    # Fix total cells
    if(!is.null(total_cells)){
	if(length(ds.cells) > total_cells){
            ds.cells = sort(resample(ds.cells, total_cells))
	}
    }

    return(ds.cells)
}


knn_downsample = function(counts, pca.rot=NULL, ident=NULL, ngene=NULL, cells.use=NULL, k1=10, k2=50, transitions=1, num_pcs=10, uniform_rates=FALSE, total_cells=NULL, cells_per_ident=NULL){

    # get input arguments
    if(is.null(ident)){ident = rep('A', ncol(counts)); names(ident) = colnames(counts)}
    if(is.null(ngene)){ngene = colSums(counts > 0); names(ngene) = colnames(counts)}

    # fix types
    if(is.null(names(ident))){names(ident) = colnames(counts)}
    if(is.null(names(ngene))){names(ngene) = colnames(counts)}
    ident = as.factor(ident)

    # subset data
    if(is.null(cells.use)){cells.use = colnames(counts)}
    counts = counts[,cells.use]
    ident = ident[cells.use]
    ngene = ngene[cells.use]

    # calculate pca
    if(is.null(pca.rot)){
        var_genes = get_var_genes(data=counts, ident=ident, method='loess', num_genes=1500)
        data = log2(calc_tpm(counts=counts) + 1)
        pca.rot = as.data.frame(rpca(t(data[var_genes,]), k=num_pcs, retx=TRUE)$x)
    }
    pca.rot = pca.rot[cells.use, 1:num_pcs]

    cat('\n\n----------\ndownsample\n----------\n')
    cat('\n\nSizes before:\n')
    print(sort(table(ident)))

    # select number of cells per group
    num_cells = sort(table(ident))
    print(table(ident))
    num_keep = num_cells_per_group(ident, total_cells, cells_per_ident)
    cat('\n\nSizes after:\n')
    print(num_keep)

    # cells to keep
    cells = data.frame(rep(NA, ncol(counts)), row.names=colnames(counts))
    rates = data.frame(rep(NA, ncol(counts)), row.names=colnames(counts))

    # downsample each group separately
    ident = droplevels(ident)
    for(group in levels(ident)){cat(paste0('\n\nDownsampling ', group, ':\n'))

        # subset data
        cells.use = colnames(counts)[ident == group]

        # skip small groups
	if(num_cells[[group]] <= num_keep[[group]]){
	    cat(paste('Selected', length(cells.use), 'cells\n'))
	    rates[cells.use, 1] = NA
	    cells[cells.use, 1] = 1
	    next
	}

	if(uniform_rates == FALSE){

	    # calculate k1 nearest neighbors
	    cat(paste('Calculating', k1, 'nearest neighbors\n'))
	    g = nng(pca.rot[cells.use,], k=k1, method='cosine', use.fnn=TRUE)

	    # calculate transition matrix
	    cat(paste('Calculating', transitions, 'transitions\n'))
	    G = t(as_adj(g))
	    G = scale(G, center=F, scale=colSums(G))
	    G = G %^% transitions

	    # number of neighbors for each cell
	    p = rowSums(G > 0) + 1

	    # solve optimal downsampling rates
	    cat('Solving optimal downsampling rates\n')
	    num_remove = num_cells[[group]] - num_keep[[group]]
	    f = function(a){(sum((exp(a)/p)^(1/(p-1))) - num_remove)^2}
	    a = optim(0, f, method='L-BFGS-B', upper=log(min(p)))$par
	    cat(paste0('lambda = ', signif(exp(a),3), ', f(lambda) = ', signif(f(a),3), '\n'))

	    d = structure((exp(a)/p)^(1/(p-1)), names=cells.use)

	} else {
	    d = (num_cells[[group]] - num_keep[[group]])/num_cells[[group]]
	    d = structure(rep(d, length(cells.use)), names=cells.use)
	}

	cat('Deletion rates:\n')
	print(quantile(d))

	# calculate k2 nearest neighbors
	cat(paste('Local subsampling with', k2, 'nearest neighbors'))
	cat('Constructing nng')
	g = nng(pca.rot[cells.use,], k=k2, method='cosine', use.fnn=TRUE)

	# downsample cells
	cat('Downsampling cells')
	keep = sapply(1:length(cells.use), function(i){
	    ci = cells.use[[i]]
	    nn = cells.use[g[[i]][[1]]]
	    return(sum(ngene[[ci]] >= ngene[nn]) >= d[[i]]*k2)
	})
	cells[cells.use,1] = keep
	rates[cells.use,1] = d
    }
    cat(paste('\n\nSelected', sum(cells[,1]), 'cells\n\n\n'))

    return(list(cells=cells, rates=rates))
}


