load_map = function(fn, key, value){

    # Construct a map between two comma-separated columns of a mapping file

    # Load mapping file
    x = read.table(fn, stringsAsFactors=F, header=T)

    # Map names to indices
    n2i = structure(1:ncol(x), names=colnames(x))

    # Construct map
    x = unique(do.call(rbind, apply(x, 1, function(xi){
        k = unlist(strsplit(xi[[n2i[[key]]]], ','))
	v = unlist(strsplit(xi[[n2i[[value]]]], ','))
	expand.grid(k,v)
    })))
    colnames(x) = c('k', 'v')

    # Aggregate by key
    x = data.frame(aggregate(v ~ k, x, function(a){c(as.character(a))}), row.names=1, stringsAsFactors=F)
    x = structure(x[,1], names=rownames(x))

    # Return a named list of keys -> values
    return(x)
}


select_keys = function(map, keys, invert=FALSE){

    # Find keys with matching or mismatching values
    # Examples:
    # > anno = load_anno(key='name', value='ident')
    #
    # > select_keys(anno, 'Monocytes')
    # - All descendent nodes:
    # - DCs, Macrophages, Monocytes (not Myeloid)
    #
    # > select_keys(anno, 'Monocytes', invert=T)
    # - Non-descendent non-ancestral nodes
    # - T cells, B cells, etc. (not Myeloid)

    i = sapply(map, function(values){
        if(invert == FALSE){
	    all(values %in% unlist(map[keys]))
	} else {
	    all(! values %in% unlist(map[keys]))
	}
    })

    return(names(map)[i])
}


flatten_keys = function(map, keys){

    # Given a map and a list of keys, this function returns:
    # keys + all non-overlapping "base" keys (length = 1)
    # Example:
    # > anno = load_anno(key='name', value='ident')
    # > good = c('Tcell', 'B', 'N', 'F_Fib', 'M')
    # > flatten_keys(anno, good)
    # > c('B', 'N', ..., 'E_Enterocytes', 'E_Tuft', ..., 'F_Glial')

    unmapped = keys[! keys %in% names(map)]
    if(length(unmapped) > 0){stop(unmapped)}
    others = select_keys(map, keys, invert=TRUE)
    others = others[sapply(map[others], length) == 1]
    return(c(keys, others))
}


basest_key = function(map, keys, multi=FALSE){

    # Find the "basest" key in a tree-structured map
    # Example:
    # anno = load_anno(key='name', value='ident')
    # keys = c('Tuft', 'Secretory', 'Epithelial')
    # basest_key(keys) = 'Tuft'

    # Test data structure
    keys = as.character(keys)
    if(length(keys) == 1){return(keys)}
    test = sapply(map[keys], function(a) sapply(map[keys], function(b) all(a %in% b)))
    diag(test) = FALSE

    # Return basest group
    if(all(test) == TRUE){
        names(which.min(sapply(map[keys], length)))
    } else if(multi == TRUE) {
        keys[! apply(test, 1, any)]
    } else {
        stop('map not hierarchical')
    }
}


highest_key = function(map, keys){

    # Find the "highest" key in a tree-structured map
    # Example:
    # anno = load_anno(key='name', value='ident')
    # keys = c('Tuft', 'Secretory', 'Epithelial')
    # highest_key(keys) = 'Epithelial'

    # Test data structure
    keys = as.character(keys)
    test = any(sapply(map[keys], function(a) sapply(map[keys], function(b) all(b %in% a))))

    # Return highest group
    if(test == TRUE){
        names(which.max(sapply(map[keys], length)))
    } else {
        stop('map not hierarchical')
    }
}


load_anno = function(type='gut', fn=NULL, key='ident', value='name', fix=c('desc', 'description', 'final')){

    # Get filename
    if(is.null(fn)){
	if(type == 'gut'){fn = '/broad/smillie-data/data/scrna/Colon-UC-Regev/anno/all.anno.txt'} else {fn = '~/ens/analysis/1.atlas/current/anno/neur.anno.txt'}
    }

    # Load annotation map
    x = load_map(fn, key, value)

    # Reformat text
    fix_text = function(a){
        a = gsub('\\\\n', '\n', a)
	a = gsub('_', ' ', a)
	a
    }
    x = sapply(x, function(a) gsub('^ ', '', a), simplify=F)

    if(key %in% fix){
        names(x) = fix_text(names(x))
    }

    if(value %in% fix){
        x = sapply(x, fix_text)
    }

    return(x)
}

