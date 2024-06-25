# source('/broad/smillie-data/code/csmillie/single_cell/map_gene.r')


psi_log = function(x, base=2, zero=NULL){
    if(is.null(zero)){
        if(any(x == 0)){zero = .5*min(x[x > 0])} else {zero = 0}
    }
    x[x == 0] = zero
    log(x, base=base)
}

geom_mean = function(x, base=2, zero=NULL){
    base**sum(psi_log(x, base=base, zero=zero))/length(x)
}

colGmeans = function(x, base=2){
    # Calculate geometric mean of gene expression for each cell
    # x = [genes x cells] matrix

    # Estimate separate zero for each gene
    zeros = apply(x, 1, function(a) if(all(a > 0)){0} else {.5*min(a[a > 0])})
    x = x + zeros
    apply(x, 2, geom_mean, base=base, zero=0)
}

# predict_dge_type = function(x, bases.use=c(2, exp(1), 10), tol=1e-4){
#     # Returns counts, tpm, or log base

#     # Select data
#     u = as.numeric(x[,1])
#     v = as.numeric(x[,2])

#     # Check for integers
#     if(abs(sum(u - round(u))) <= tol & abs(sum(v - round(v))) <= tol){
#         return('counts')
#     }

#     # Check for constant sum
#     if(abs(sum(u) - sum(v)) <= tol){
#         return('tpm')
#     }

#     # Check for big numbers
#     if(max(u) > 32 | max(v) > 32){
#         print('predict_log_base: found big numbers, so guessing data type = imputed TPM')
#         return('tpm')
#     }

#     # Check for logTPM
#     u = sapply(bases.use, function(base){
#         abs(sum(base**u) - sum(base**v))
#     })
#     i = which.min(u)

#     # Test logTPM deviance
#     if(u[[i]] <= tol){
#         return(bases.use[[i]])
#     } else {
#         print('predict_log_base: confused, so guessing data type = TPM???')
#     }
# }


calc_tpm = function(obj=NULL, counts=NULL, tpm=NULL, data=NULL, genes.use=NULL, cells.use=NULL, total=1e4){
    # Calculate TPM from $counts or data argument

    # Select counts data
    if(!is.null(counts)){data = counts; type = 'counts'}
    if(!is.null(tpm)){data = tpm; type = 'tpm'}
    if(is.null(data)){
        print('Guessing data type = counts')
        data = obj$counts
	type = 'counts'
    }
    if(is.null(type)){
        type = predict_dge_type(data[, sample(1:ncol(data), 2)])
    }
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}

    # Get genes and cells
    genes.use = intersect(genes.use, rownames(data))
    cells.use = intersect(cells.use, colnames(data))

    if(type == 'counts'){
        require(wordspace)
        numi = colSums(data[,cells.use])
	data = scaleMargins(data[genes.use, cells.use, drop=F], cols=total/numi)
    }
    if(type == 'tpm'){
        # Do nothing
    }
    if(is.numeric(type)){
        data = type**data[genes.use, cells.use, drop=F]
	data = data - min(data)
    }

    return(data)
}


# get_data = function(obj, data.use='tpm', tpm=NULL, genes.use=NULL, cells.use=NULL, qnorm=FALSE){

#     # Retrieve data from a singlecell object
#     # data.use can be: counts, tpm, log2, data, any matrix
#     # optionally, pre-calculate tpm for speed

#     if(!is.character(data.use)){
#         data = data.use
# 	qnorm = FALSE
#     } else if(data.use == 'counts'){
#         data = obj$counts
#     } else if(data.use == 'tpm'){
#         if(is.null(tpm)){
# 	    data = calc_tpm(obj, genes.use=genes.use, cells.use=cells.use)
# 	} else {
#             data = tpm
#         }
#     } else if(data.use == 'log2'){
#         if(predict_dge_type(obj$data, bases.use=c(2)) != 2){
# 	    print('Warning: obj$data is not log base 2')
# 	}
# 	data = obj$data
#     } else if(data.use == 'scale'){
#         if(!is.null(genes.use)){
# 	    genes.use = intersect(genes.use, rownames(obj$data))
# 	} else {
# 	    genes.use = rownames(obj$data)
# 	}
# 	data = t(scale(t(obj$data[genes.use,])))
#     } else {
#         stop('Error: get_data invalid argument')
#     }

#     # Quantile normalization
#     if(qnorm == TRUE){
#         data = quantile_normalize(data=data)
#     }
#     data = as(data, 'sparseMatrix')

#     # Subset genes and cells
#     if(!is.null(genes.use)){data = data[intersect(genes.use, rownames(data)),]}
#     if(!is.null(cells.use)){data = data[,intersect(cells.use, colnames(data))]}

#     return(data)
# }


# map_names = function(obj=NULL, data=NULL, meta=NULL, names=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, source='auto', target='auto'){

#     # Map input arguments to genes and feats

#     # Initialize variables
#     names = as.list(names)
#     genes = feats = c()
#     regex = as.list(regex)

#     # Get data and metadata
#     if(is.null(data)){data = obj$data}
#     if(is.null(meta)){meta = obj$meta.data}

#     # Map names
#     if(!is.null(names)){
#         if(target == 'auto'){target = predict_organism(rownames(data)[1:100])}
# 	genes = sapply(names, function(a){a = as.character(a)
# 	    i = a %in% rownames(data)
# 	    u = a[i]
# 	    if(sum(!i) > 0){u = c(u, intersect(map_gene(a[!i], source=source, target=target), rownames(data)))}
# 	    unique(u)
# 	}, simplify=F)
# 	feats = sapply(names, function(a) a[a %in% colnames(meta)])
#     }

#     # Map regex
#     if(!is.null(regex)){ # wrap in structure? names = regex
#         genes = c(genes, sapply(regex, function(a) grep(a, rownames(data), value=T, perl=T), simplify=F))
# 	feats = c(feats, sapply(regex, function(a) grep(a, colnames(meta), value=T, perl=T), simplify=F))
#     }

#     # Map files
#     if(!is.null(files)){
#         sig = do.call(c, lapply(files, function(file) load_signature(file, file.regex=file.regex, file.cols=file.cols)))
# 	sig = sapply(sig, function(a) intersect(map_gene(a, source=source, target=target), rownames(data)), simplify=F)
# 	genes = c(genes, sig)
#     }

#     # Filter genes and feats
#     genes = genes[sapply(genes, length) > 0]
#     feats = feats[sapply(feats, length) > 0]
#     if(!is.null(top)){
#         genes = sapply(genes, function(a) as.character(na.omit(a[1:top])), simplify=F)
# 	feats = sapply(feats, function(a) as.character(na.omit(a[1:top])), simplify=F)
#     }

#     # Fix gene names
#     genes = genes[sapply(genes, length) > 0]
#     if(length(genes) > 0){
#         if(is.null(names(genes))){names(genes) = rep('', length(genes))}
# 	names(genes)[names(genes) == ''] = sapply(genes[names(genes) == ''], paste, collapse='.')
#     }

#     # Fix feat names
#     feats = feats[sapply(feats, length) > 0]
#     if(length(feats) > 0){
#         if(is.null(names(feats))){names(feats) = rep('', length(feats))}
# 	names(feats)[names(feats) == ''] = sapply(feats[names(feats) == ''], paste, collapse='.')
#     }

#     return(list(genes=genes, feats=feats))
# }


# score_cells = function(obj=NULL, names=NULL, data=NULL, meta=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, source='auto', target='auto', scores=NULL,
#                        data.use='log2', combine_genes='mean', groups=NULL, group_stat='mean', genes_first=TRUE, cells.use=NULL, make.names=TRUE, do.log=FALSE, drop_zeros=TRUE, qnorm=FALSE){

#     require(Matrix)
#     require(Matrix.utils)

#     # Score gene expression across cells and optionally aggregate
#     # The default steps are:
#     # 1) Select data (data.use = 'tpm', 'log2', 'scale')
#     # 2) Calculate mean expression across genes (combine_genes = 'sum', 'mean', 'geom', 'max', 'norm', 'scale', 'scale2')
#     # 3) Calculate mean expression within each cell type (group_by, group_stat = 'mean', 'alpha', 'mu')
#     # 4) Log transform results (do.log)
#     # If genes_first == FALSE, then calculate the group means *before* combining across genes

#     # Fix input arguments and get data for name mapping
#     if(is.null(data)){map_data = obj$data} else {map_data = data}
#     if(is.null(meta)){meta = obj$meta.data}
#     if(!is.null(groups)){names(groups) = colnames(obj$data)}
#     if(!is.null(scores)){
#         scores = as.data.frame(scores, stringsAsFactors=F)
# 	if(length(intersect(rownames(scores), colnames(obj$data))) == 0){
# 	    if(nrow(scores) == length(cells.use)){
# 	        rownames(scores) = cells.use
# 	    } else if(nrow(scores) == ncol(obj$data)){
# 	        rownames(scores) = colnames(obj$data)
# 	    } else {
# 	        stop('Cannot find rownames for scores')
# 	    }
# 	}
# 	scores = scores[colnames(obj$data),,drop=F]
#     }

#     # Map genes and feats
#     res = map_names(obj=obj, data=map_data, meta=meta, names=names, regex=regex, files=files, file.cols=file.cols, file.regex=file.regex, top=top, source=source, target=target)
#     genes = res$genes
#     feats = res$feats
#     genes.use = unique(do.call(c, genes))
#     if(length(genes) == 0 & length(feats) == 0){return(scores)}

#     # Select data with genes.use
#     if(is.null(data)){data = get_data(obj, data.use=data.use, cells.use=cells.use, qnorm=qnorm)}

#     # Subset cells
#     if(!is.null(cells.use)){
#         data = data[,cells.use,drop=F]
#         meta = meta[cells.use,,drop=F]
#         if(!is.null(groups)){groups = groups[cells.use]}
# 	if(!is.null(scores)){scores = scores[cells.use,,drop=F]}
#     }

#     group_genes = function(x, method){

#         # combine expression data across genes within a signature
# 	# x = [genes x cells] matrix
# 	# method = 'sum', 'mean', or 'gmean'
# 	# returns [genes x cells] or [1 x cells] matrix

# 	if(nrow(x) == 1){return(x[1,,drop=F])}
# 	if(method == 'sum'){
# 	    t(colSums(x))
# 	} else if(method == 'mean'){
# 	    t(colMeans(x, na.rm=T))
# 	} else if(method == 'geom'){
# 	    t(colGmeans(x))
# 	} else if(method == 'max'){
# 	    x = x/apply(x, 1, max)
# 	    t(colMeans(x, na.rm=T))
# 	} else if(method == 'norm'){
# 	    x = x/apply(x, 1, sum)
# 	    t(colMeans(x, na.rm=T))
# 	} else if(method == 'scale'){
# 	    x = t(scale(t(x)))
# 	    t(colMeans(x, na.rm=T))
# 	} else if(method == 'scale2'){
# 	    x = t(scale(t(x), center=F))
# 	    t(colMeans(x, na.rm=T))
# 	} else if(method == 'none'){
# 	    x
# 	}else {
# 	    stop('Error: invalid combine_genes method')
# 	}
#     }

#     group_cells = function(x, groups, method){

#         # combine expression data across cells
# 	# x = [genes x cells] matrix
# 	# group_stat = 'alpha', 'mu', or 'mean'
# 	# returns [genes x groups] matrix

# 	if(is.null(groups)){return(x)}
# 	if(method %in% c('n', 'sum')){
# 	    if(method == 'n'){x = x > 0}
# 	    x = t(data.frame(aggregate(t(x), list(groups), sum, na.rm=T), row.names=1))
# 	} else {
# 	    if(method == 'alpha'){x = x > 0}
# 	    if(method == 'mu'){x[x == 0] = NA}
# 	    x = t(data.frame(aggregate(t(x), list(groups), mean, na.rm=T), row.names=1))
# 	}
# 	x[is.na(x)] = 0
# 	x
#     }

#     # Calculate scores
#     names.use = unique(c(names(genes), names(feats)))

#     # Fast indexing for flat structures
#     name_map = sapply(names.use, function(a) c(genes[[a]], feats[[a]]), simplify=F)
#     do.flat = all(sapply(name_map, length) == 1)
#     if(do.flat == TRUE){
#         genes[['flat']] = do.call(c, genes)
# 	feats[['flat']] = do.call(c, feats)
# 	names.iter = 'flat'
# 	combine_genes = 'none'
#     } else {
#         names.iter = names.use
#     }

#     backup = scores
#     scores = lapply(names.iter, function(name){

#         # Combine data and metadata
# 	if(name %in% names(genes)){si = data[genes[[name]],,drop=F]} else {si = c()}
# 	if(name %in% names(feats)){if(is.null(si)){si = t(meta[,feats[[name]],drop=F])} else {si = rBind(si, t(meta[,feats[[name]],drop=F]))}}
# 	si = as.matrix(si)

# 	if(genes_first == TRUE){
# 	    si = group_genes(si, method=combine_genes)
# 	    si = group_cells(si, groups=groups, method=group_stat)
# 	} else {
# 	    si = group_cells(si, groups=groups, method=group_stat)
# 	    si = group_genes(si, method=combine_genes)
# 	}

# 	if(do.log == TRUE){
# 	    si = psi_log(si, base=2)
# 	}
# 	si = data.frame(t(si), stringsAsFactors=F)
#     })

#     # Collapse scores
#     if(do.flat == TRUE){
# 	scores = scores[[1]][,make.names(name_map[names.use]),drop=F]
#     } else {
#         do.collapse = all(lapply(scores, ncol) == 1)
# 	if(do.collapse == TRUE){
# 	    scores = as.data.frame(do.call(cbind, scores), stringsAsFactors=F)
# 	}
#     }

#     # Fix names
#     if(make.names == TRUE){names.use = make.names(names.use)}
#     names(scores) = names.use
#     if(drop_zeros == FALSE){cols.use = names.use; if(is.null(cols.use)){cols.use=names.use}; scores[,setdiff(cols.use, colnames(scores))] = 0; scores = scores[,cols.use]}

#     # Combine data
#     if(!is.null(backup)){
#         if(is.data.frame(scores)){
#             scores = cbind.data.frame(scores, backup)
#         } else {
#             scores = c(scores, backup)
#         }
#     }

#     return(scores)
# }


# score_signatures = function(obj=NULL, gene_sets=NULL, cells.use=NULL, cells.bg=NULL, data.use='log2', g=20, nperm=100, do.collapse=TRUE){
#     library(Hmisc)

#     # fix input arguments
#     data = get_data(obj, data.use=data.use, cells.use=cells.use)
#     if(!is.list(gene_sets)){gene_sets = list(gene_sets)}
#     if(is.null(cells.bg)){
#         cells.bg = colnames(data)
#     } else {
#         cells.bg = intersect(cells.bg, colnames(data))
#     }

#     # divide genes into g equal-frequency expression bins
#     g2b = setNames(cut2(rowMeans(data[,cells.bg]), g=g), rownames(data))
#     b2g = sapply(levels(g2b), function(a) names(g2b)[g2b == a])

#     # iterate over gene sets
#     sapply(names(gene_sets), function(name.use){print(name.use)

#         tryCatch({

#         genes.use = intersect(gene_sets[[name.use]], rownames(data))

#     	# get background gene sets
#     	genes.bg = sapply(genes.use, function(gene){
#             bin = g2b[[gene]]
# 	    sample(b2g[[bin]], nperm, replace=TRUE)
#         })

#         # calculate signature scores
#         true = t(data[genes.use,])
#         null = apply(genes.bg, 2, function(a) colMeans(data[a,]))

#         # average score across genes
#     	if(do.collapse == TRUE){
#             true = rowMeans(true)
# 	    null = rowMeans(null)
# 	    true - null
#         } else {
#             list(true=true, null=null)
#         }

# 	}, error = function(e){
# 	    setNames(rep(NA, ncol(data)), colnames(data))
# 	})

#     }, simplify=F)
# }

nice_agg = function(x, g, type='median'){
    # memory efficient version of data.frame(aggregate(x, list(g), mean), row.names=1)
    if(nrow(x) != length(g)){
        stop('error: invalid dimensions (transpose data?)')
    }
    g = as.factor(g)
    y = sapply(levels(g), function(gi){
        i = (g == gi)
	i[is.na(i)] = FALSE
	if(sum(i) > 0){
	    if(type == 'mean'){
	        colMeans(x[i,,drop=F])
	    } else if(type == 'sum'){
	        colSums(x[i,,drop=F])
	    } else if(type == 'sd'){
	        apply(x[i,,drop=F], 2, sd)
	    } else if(type == 'median'){
	        apply(x[i,,drop=F], 2, median)
	    }
	} else {
	    rep(NA, ncol(x))
	}
    })
    y = t(y)
    rownames(y) = levels(g)
    colnames(y) = colnames(x)
    y
}
