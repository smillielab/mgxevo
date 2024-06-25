run_dmap = function(data, n_pcs=0, dpt=FALSE, root_cell=NULL, n_neighbors=10, paga=FALSE, clusters=NULL, cleanup=TRUE, out=NULL){

    # -----------------------------------------
    # Run scanpy on [cells x feats] data matrix
    # -----------------------------------------

    # Write data
    if(is.null(out)){
        out = tempfile(pattern='scanpy.', tmpdir='~/tmp', fileext='.data.txt')
    } else {
        out = paste0(out, '.data.txt')
    }
    write.table(data, file=out, sep='\t', quote=F, row.names=FALSE, col.names=FALSE)

    # Write clusters
    if(is.null(clusters)){clusters = 'louvain_groups'} else {
        clusters_fn = gsub('data', 'clusters', out)
	write.table(clusters, clusters_fn, quote=F, row.names=FALSE, col.names=FALSE)
    }

    # Root cell index
    if(is.null(root_cell)){iroot = 0} else {iroot = match(root_cell, rownames(data))-1}

    # Diffusion map
    command = paste('python ~/code/single_cell/run_dmap.py --data', out, '--n_pcs', n_pcs, '--n_neighbors', n_neighbors, '--iroot', iroot, '--out', gsub('.data.*', '', out))

    # Pseuodotime
    if(dpt == TRUE){command = paste(command, '--dpt')}

    # Approximate graph abstraction
    if(paga == TRUE){command = paste(command, '--paga', '--clusters', clusters_fn)}

    system(command)

    # Load & cleanup results
    res = list(dmap='X_diffmap', dpt_time='dpt_pseudotime', paga_full='connectivities', paga_tree='connectivities_tree', categories='categories')

    for(name in names(res)){
        fn = gsub('data', res[[name]], out)
	if(file.exists(fn)){
	    res[[name]] = as.data.frame(fread(fn, header=FALSE))
	    if(nrow(data) == nrow(res[[name]])){
	        rownames(res[[name]]) = rownames(data)
	    }
	    if(cleanup == TRUE){system(paste('rm', fn))}
	}
    }
    tryCatch({for(name in c('paga_tree', 'paga_full')){rownames(res[[name]]) = colnames(res[[name]]) = res[['categories']][,1]}}, error=function(e){})
    if(cleanup == TRUE){system(paste('rm', out))}

    return(res)
}
