select_bg_genes = function(obj=NULL, scores=NULL, genes.use, cells.use=NULL, nbins=20, nperm=100, type='mean', cutoff=0.25){
    library(Hmisc)

    # ---------------------------------------------------------------------------------------------------
    # Given singlecell object and target gene set, select [nperm] expression-matched background gene sets
    # ---------------------------------------------------------------------------------------------------
    # type:
    # - mean = match mean expression level across all cells
    # - num_ident = match number of idents with minimum expression cutoff (e.g. number of cell types with mean expression > 0.25)

    # fix input arguments
    if(is.null(cells.use)){
        cells.use = colnames(obj$data)
    }

    # divide genes into g equal-frequency expression bins
    if(is.null(scores)){
        if(type == 'mean'){
            scores = rowMeans(obj$data[,cells.use])
        } else {
            scores = nice_agg(t(obj$data), obj$ident)
	    scores = colSums(scores >= cutoff)
        }
    }
    genes2bin = setNames(cut2(scores, g=nbins), names(scores))
    bin2genes = sapply(levels(genes2bin), function(a) names(genes2bin)[genes2bin == a])

    # select control genes
    genes.use = intersect(genes.use, names(scores))

    # background gene sets
    genes.bg = sapply(genes.use, function(gene){
        bin = genes2bin[[gene]]
	sample(bin2genes[[bin]], nperm, replace=TRUE)
    })

    if(nperm > 1){genes.bg = t(genes.bg)}
    return(genes.bg)
}

tm_modules = function(obj, genes.use=NULL, cells.use=NULL, ident.use=NULL, cells_per_group=NULL, binarize=FALSE){

    # fix input arguments
    if(is.null(genes.use)){
        genes.use = rownames(obj$data)
    }
    if(is.null(cells.use)){
        cells.use = colnames(obj$data)
    }
    if(is.null(ident.use)){
        ident.use = obj$ident
    }
    if(!is.null(cells_per_ident)){
        cells.use = simple_downsample(colnames(obj$data), ident.use, cells_per_group=cells_per_group)
    }

    # construct document-term-matrix
    data.use = obj$counts[genes.use, cells.use]
    data.use = as(data.use, 'dgTMatrix')
    data.use = as.simple_triplet_matrix(data.use)

    # binarize frequencies (optional)
    if(binarize == TRUE){
        dtm$Frequency = 1
    }

}


get_module = function(obj, gene){
    u = as.numeric(obj$data[gene,])
    v = t(as.matrix(obj$data))
    sort(cor(u, v)[1,])
}

gene_modules = function(obj, ident, data, gene, num_pcs=10, method='pearson', scale=FALSE){

    require(rsvd)
    require(ppcor)

    # subsample data
    cells.use = intersect(rownames(data), colnames(obj$data)[obj$ident == ident])
    ngene = scale(obj$meta.data[cells.use, 'nGene'])

    # memory efficient scaling
    data = data[cells.use,]
    if(scale == TRUE){
        g = split(1:ncol(data), cut(1:ncol(data), 10))
	for(j in g){
	    data[,j] = scale(data[,j])
	}
    }
    data = as.data.frame(data)

    # type 1 = correlation
    mod1 = sort(cor(data[,gene], data, method=method)[1,])

    # type 2 = partial correlation
    mod2 = sort(sapply(data, function(a) tryCatch(pcor.test(data[,gene], a, ngene, method=method)$estimate, error=function(e) NA)))

    # calculate pca
    pca.obj = rpca(data, center=FALSE, scale=FALSE, retx=TRUE, k=20)

    # type 3 = gene loadings
    mod3 = t(pca.obj$rotation)
    mod3 = sort(cor(mod3[,gene], mod3, method=method)[1,])

    # set up linear model
    x = as.data.frame(pca.obj$x[,1:num_pcs])
    colnames(x) = paste0('PC', 1:ncol(x))
    x$nGene = ngene
    x = as.matrix(x)

    # type 4 = correlation of residuals
    r = sapply(data, function(a) .lm.fit(x, a)$residuals)
    mod4 = sort(cor(r[,gene], r, method=method)[1,])

    # type 5 = linear models
    mod5 = sort(sapply(data, function(a) .lm.fit(cbind(a, x), data[,gene])$coefficients[[1]]))

    list(mod1, mod2, mod3, mod4, mod5)
}


mi_modules = function(data, genes.use, cells.use=NULL, nbins=8){
    library(data.table)
    library(entropy)

    # fix data
    data = as.data.frame(data)
    if(!is.null(cells.use)){
        cells.use = intersect(cells.use, rownames(data))
        data = data[cells.use,]
    }
    if(nrow(data) == 0){
        stop('error: cell names')
    }
    if(sum(is.na(data)) > 0){
        stop('error: missing values')
    }

    # discretize
    data = data[,sapply(data, sd) > 0]
    data = lapply(data, function(a) cut(a, 8, include.lowest=T))

    # pairwise mutual information
    sapply(genes.use, function(a) sapply(data, function(b) tryCatch({mi.empirical(table(data[[a]], b))}, error=function(e){NA})))
}
