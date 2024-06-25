require(Matrix)
require(methods)
require(data.table)
require(tibble)
require(naturalsort)
require(wordspace)

options(max.print=1000)

source('batch.r')
source('cluster.r')
source('contamination.r')
source('dmap.r')
source('downsample.r')
source('frequencies.r')
source('gsea.r')
source('map.r')
source('markers.r')
source('parallel.r')
source('pca.r')
source('plot.r')
source('specific.r')
source('tsne.r')
source('var_genes.r')
source('../../ut/mtx.r')
# source('/broad/smillie-data/code/csmillie/FIt-SNE/fast_tsne.R')


msg = function(text, verbose){
    if(verbose == TRUE){
        print(text)
    }
}


set.ident = function(obj, ident.use){
    obj$ident = setNames(as.factor(ident.use), colnames(obj$data))
    return(obj)
}


autoload_obj = function(prefix, ming=0, minc=0){

    # read data
    data = paste0(prefix, '.matrix.mtx')
    if(file.exists(data)){
        data = readMM(data)
    } else {
        data = paste0(prefix, '.matrix.csv')
	data = read.table(data, sep=',')
    }
    print(paste('data', dim(data)))

    # read genes (rownames)
    rows = paste0(prefix, '.genes.csv')
    if(file.exists(rows)){
        rows = read.table(rows, sep=',')[,1]
    }

    # read barcodes (colnames)
    cols = paste0(prefix, '.barcodes.csv')
    if(file.exists(cols)){
        cols = read.table(cols, sep=',')[,1]
    }

    # fix rownames/colnames
    if(abs(nrow(data) - length(cols)) < abs(nrow(data) - length(rows))){
        data = t(data)
    }
    if(nrow(data) == (length(rows) - 1)){
        rows = rows[2:length(rows)]
    }
    if(ncol(data) == (length(cols) - 1)){
        cols = cols[2:length(cols)]
    }

    # set rownames/colnames
    rownames(data) = rows
    colnames(data) = cols

    # read metadata
    meta = paste0(prefix, '.meta.csv')
    if(file.exists(meta)){
        meta = read.table(meta, sep=',', header=TRUE, row.names=1)
	print(paste('metadata:', dim(meta)))
    }

    # read image
    image = paste0(prefix, '.image_hires.png')
    if(file.exists(image)){
	library(png)
        image = readPNG(image)
	print(paste('image:', dim(image)))
    }

    # make singlecell object
    obj = make_obj(counts=data, minc=minc, ming=ming)
    obj$meta = meta
    obj$image = image

    return(obj)
}


init_obj = function(){
    obj = list()
    obj$counts = NA
    obj$data = NA
    obj$meta.data = NA
    obj$ident = NA
    obj$tsne.rot = NA
    obj$image = NA
    obj$pca.obj = NA
    obj$pca.rot = NA
    return(obj)
}

seur2obj = function(seur){

    obj = init_obj()
    
    if('raw.data' %in% slotNames(seur)){
        obj$counts = seur@raw.data; seur@raw.data = data.frame()
        obj$data = seur@data; seur@data = data.frame()
    	obj$meta.data = seur@data.info; seur@data.info = data.frame()
    	obj$ident = seur@ident; seur@ident = NA
    	obj$tsne.rot = seur@tsne.rot; seur@tsne.rot = data.frame()
    	obj$pca.obj = seur@pca.obj; seur@pca.obj = list()
    	obj$pca.rot = seur@pca.rot; seur@pca.rot = data.frame()
    }
    
    if('assays' %in% slotNames(seur)){
        obj$counts = seur@assays$RNA@counts
    	obj$data = seur@assays$RNA@data
    	obj$meta.data = seur@meta.data
    	obj$ident = seur@active.ident
    	obj$tsne.rot = as.data.frame(seur@reductions$umap@cell.embeddings); colnames(obj$tsne.rot) = c('tSNE_1', 'tSNE_2')
    	obj$pca.obj = seur@reductions$pca
    	obj$pca.rot = seur@reductions$pca@cell.embeddings
    }

    seur = data.frame()
    return(obj)
}

make_obj = function(counts=NULL, regex='', regexv='', minc=10, maxc=NULL,maxc_per_group=NULL, ming=500, maxg=1e6, genes.use=NULL, cells.use=NULL, ident_fxn=NULL, verbose=FALSE, x11=FALSE, bayes.n=0, qnorm=F){

    # Load packages
    require(data.table)
    require(Matrix)
    source('/broad/smillie-data/code/csmillie/util/mtx.r')

    # Set graphics device
    options(device=pdf)

    # Load counts from singlecell object, matrix, or file
    msg('loading counts', verbose)
    if(typeof(counts) == typeof('')){
        if(file.exists(counts)){
	    counts = fread(paste('zcat', counts))
	    counts = data.frame(counts, row.names=1)
	} else {
	    counts = read_mtx(prefix=counts)
	}
    }
    msg( sprintf('Counts = %d x %d', nrow(counts), ncol(counts)), verbose)

    # Subset counts with NA, regex, genes.use, and cells.use
    msg( 'Subsetting counts', verbose)
    if(regex != ''){
        j = grep(regex, colnames(counts))
	counts = counts[,j,drop=F]
    }
    if(regexv != ''){
        j = grep(regexv, colnames(counts), invert=T)
	counts = counts[,j,drop=F]
    }
    if(!is.null(genes.use)){
	genes.use = intersect(rownames(counts), genes.use)
	counts = counts[genes.use,,drop=F]
    }
    if(!is.null(cells.use)){
	cells.use = intersect(colnames(counts), cells.use)
	counts = counts[,cells.use,drop=F]
    }
    genes.use = rowSums(is.na(counts)) == 0
    counts = counts[genes.use,,drop=F]

    msg( sprintf('counts = %d x %d', nrow(counts), ncol(counts)), verbose)
    
    # Convert counts to sparse matrix
    if(is.data.frame(counts)){counts = as.matrix(counts)}
    counts = as(counts, 'sparseMatrix')
    
    # Get cell identities
    if(is.null(ident_fxn)){
        ident = sapply(strsplit(colnames(counts), '\\.'), '[[', 1)
    } else {
        ident = sapply(colnames(counts), ident_fxn)
    }
    ident = setNames(ident, colnames(counts))

    # Downsample cells
    if(!is.null(maxc) | !is.null(maxc_per_group)){
        cells.use = simple_downsample(cells=colnames(counts), groups=ident, total_cells=maxc, cells_per_group=maxc_per_group)
	msg( paste('Downsampling to', length(cells.use), 'total cells'), verbose)
	counts = counts[,cells.use]
    }

    # Filter cells by minc, ming, maxg
    msg( 'Filtering counts', verbose)
    j1 = colSums(counts > 0) >= ming
    j2 = colSums(counts > 0) <= maxg
    counts = counts[,(j1 & j2)]
    i = rowSums(counts > 0) >= minc
    counts = counts[i,]
    msg( sprintf('counts = %d x %d', nrow(counts), ncol(counts)), verbose)

    # Add Bayesian prior
    if(bayes.n > 0){
        print('Bayesian prior')
        p = scaleMargins(counts, cols=1/colSums(counts))
	p = rowMeans(p)
	p = rmultinom(n=ncol(counts), size=bayes.n, prob=p)
	print(paste('Adding', mean(colSums(p)), 'reads to every cell'))
	counts = counts + p
    }
    
    # Make singlecell object
    msg( 'Making singlecell object', verbose)
    obj = init_obj()
    obj$counts = counts
    obj$ident = ident[colnames(counts)]
    obj$meta.data = data.frame(nGene=colSums(counts > 0), nUMI=colSums(counts))
    rownames(obj$meta.data) = colnames(counts)
    
    # Normalize data (default = TPM)
    msg( 'Normalizing data', verbose)
    if(qnorm == FALSE){
        obj$data = calc_tpm(counts=obj$counts)
	obj$data@x = log2(obj$data@x + 1)
    } else {
        obj = quantile_normalize(obj=obj)
    }
    
    # Get cell identities
    msg( 'Setting cell identities', verbose)
    if(!is.null(ident_fxn)){
        ident = sapply(colnames(obj$data), ident_fxn)
	obj = set.ident(obj, ident.use=ident)
	obj$meta.data$orig.ident = obj$ident
    }

    if(length(unique(obj$ident)) > 100){
        msg( 'WARNING: nlevels(obj$ident) > 100', verbose)
        #obj$ident = '1'
    	#obj$meta.data$orig.ident = '1'
    }

    print(table(obj$ident))
    return(obj)
}

quantile_normalize = function(obj=NULL, data=NULL, zero_cut=1e-4){
    print(paste0('Quantile normalization with zero_cut = ', zero_cut))
    library(preprocessCore)
    if(!is.null(obj)){
        print(dim(obj$counts))
        obj$counts[] = as(normalize.quantiles.robust(as.matrix(obj$counts)), 'sparseMatrix')
	obj$data[] = as(normalize.quantiles.robust(as.matrix(obj$data)), 'sparseMatrix')
        obj
    } else {
        print(dim(data))
        data[] = normalize.quantiles.robust(as.matrix(data))
	print(paste('Flooring', sum(data < zero_cut), 'zeros'))
	data[data < zero_cut] = 0
	as(data, 'sparseMatrix')
    }
}


run_seurat = function(name, obj=NULL, counts=NULL, regex='', regexv='', cells.use=NULL, genes.use=NULL, minc=5, maxc=1e6, ming=200, maxg=1e6, ident_fxn=NULL, varmet='loess', var_regexv=NULL,
             var_remove=NULL, min_cv2=.25, var_genes=NULL, use_var=FALSE, bayes.n=0, qnorm=F, num_genes=1500, regress=NULL, do.batch='none', batch.use=NULL, store.batch=FALSE, bc.data=NULL,
	     design=NULL, pc.data=NULL, num_pcs=0, pcs.use=NULL, pcs.rmv=NULL, robust_pca=F, scale.max=10, 
	     perplexity=25, max_iter=1000, dist.use='euclidean', do.largevis=FALSE, do.umap=FALSE, largevis.k=50, do.fitsne=FALSE, fitsne.K=-1,
	     cluster='infomap', k=c(), verbose=T, write_out=T, do.backup=F, ncores=1, stop_cells=50, marker.test=''){

    # check input arguments
    if(! do.batch %in% c('none', 'combat', 'mnn', 'cca', 'multicca', 'liger')){stop('do.batch must be none, combat, mnn, or cca')}
    if(! is.null(batch.use)){if(is.null(names(batch.use))){stop('batch.use needs names')}}

    # Make singlecell object
    if(is.null(obj)){
        obj = make_obj(counts=counts, regex=regex, regexv=regexv, minc=minc, maxc=maxc, ming=ming, maxg=maxg, genes.use=genes.use, cells.use=cells.use, ident_fxn=ident_fxn, verbose=verbose, qnorm=qnorm, bayes.n=bayes.n)
    }
    if(ncol(obj$data) <= stop_cells){return(obj)}
    
    msg( 'Selecting variable genes', verbose)
    ident = obj$ident
    if(is.null(var_genes)){
        gi = rownames(obj$counts)
	print(paste('Starting with', length(gi), 'genes'))
	if(!is.null(var_regexv)){gi = grep(var_regexv, gi, invert=T, value=T)}
	print(paste('var_regexv:', length(gi), 'genes'))
	if(!is.null(var_remove)){gi = setdiff(gi, var_remove)}
	print(paste('var_remove:', length(gi), 'genes'))
	var_genes = get_var_genes(obj$counts, ident=ident, method=varmet, genes.use=genes.use, num_genes=num_genes, min_ident=25, use_var=use_var)
    }
    #if(is.null(var_genes)){var_genes = get_var_genes(obj$counts, ident=ident, method=varmet, num_genes=num_genes, min_ident=25)}
    #if(!is.null(var_regexv)){var_genes = grep(var_regexv, var_genes, invert=T, value=T)}
    msg( sprintf('Found %d variable genes', length(var_genes)), verbose)
    obj$var.genes = intersect(var_genes, rownames(obj$data))
    print(var_genes)
    
    # Regression
    if(!is.null(regress)){
        print('Regressing out gene signature from pc.data')
	print(obj$data[1:5,1:5])
        x = score_cells(obj, names=regress)
	obj$data = t(apply(obj$data, 1, function(a) .lm.fit(as.matrix(x), as.matrix(a))$residuals[,1]))
	print(obj$data[1:5,1:5])
    }
    
    # Batch correction with variable genes
    if(do.batch != 'none'){
	msg( 'Batch correction', verbose)
	if(is.null(batch.use)){
	    batch.use = obj$ident
	}
	batch.use = batch.use[names(obj$ident)]
	print(table(batch.use))
	if(!is.null(design)){
	    design = design[names(obj$ident),,drop=F]
	}
	if(is.null(bc.data)){
	    bc.data = batch_correct(obj, batch.use, design=design, method=do.batch, genes.use=obj$var.genes, ndim=num_pcs)
	}
	
	# store batch corrected data
	if(store.batch == TRUE){
	    obj$bc.data = bc.data
	}
	
	# write batch corrected data to file
	if(write_out == TRUE){fwrite(as.data.table(bc.data), file=paste0(name, '.bc.data.txt'), sep='\t')}
	
	pc.data = t(scale(t(bc.data), center=F))
	#pc.data = bc.data
    }
    if(is.null(pc.data)){pc.data = obj$data}

    # Number of significant PCs (stored in obj$meta.data$num_pcs)
    if(num_pcs == 0){
	num_pcs = sig.pcs.perm(scale(t(obj$data)), randomized=T, n.cores=ncores)$r + 2
    }
    if(is.na(num_pcs)){num_pcs = 5}
    obj$meta.data$num_pcs = num_pcs
    msg( sprintf('Found %d significant PCs', num_pcs), verbose)
    
    # Fast PCA on data
    if(do.batch %in% c('liger', 'multicca')){
        obj$pca.rot = as.data.frame(pc.data) # liger and multi-cca output saved in pc.data
	print(dim(obj$pca.rot))
    } else {
        obj = run_rpca(obj, data=pc.data, k=50, genes.use=obj$var.genes, robust=robust_pca, rescale=T, scale.max=scale.max)
	if(write_out == TRUE){saveRDS(obj$pca.obj, file=paste0(name, '.pca.rds'))}
        msg( 'PC loadings', verbose)
        loaded_genes = get.loaded.genes(obj$pca.obj[[1]], components=1:num_pcs, n_genes=20)
        print(loaded_genes)
    }

    # Fix problem with duplicates
    obj$pca.rot[,num_pcs] = obj$pca.rot[,num_pcs] + runif(nrow(obj$pca.rot), min=-1e-8, max=1e-8)

    # Regress out PCs
    if(is.null(pcs.use)){
        pcs.use = 1:(min(ncol(obj$pca.rot), num_pcs))
    }

    # TSNE
    knn = NULL
    if(max_iter > 0){
    if(do.fitsne == TRUE){
        msg( 'FIt-SNE', verbose)
	q = fftRtsne(obj$pca.rot[,pcs.use], max_iter=max_iter, perplexity=perplexity, K=fitsne.K, fast_tsne_path='/broad/smillie-data/code/csmillie/FIt-SNE/bin/fast_tsne')
    } else if(do.largevis == TRUE){
        msg( 'largeVis', verbose)
	q = run_largevis(obj$pca.rot[,pcs.use], k=largevis.k, save_knn=TRUE, save_weights=FALSE, dist.use=dist.use, verbose=T)
	knn = q$knn
	q = q$coords
    } else if(do.umap == TRUE){
        msg( 'umap', verbose)
	library(umap)
	q = umap(obj$pca.rot[,pcs.use], method='umap-learn')$layout
	print(dim(q))
	#q = umap(t(as.matrix(pc.data)), method='umap-learn')$layout
    } else {
        msg( 'TSNE', verbose)
	require(Rtsne)
        if(dist.use == 'euclidean'){
            q = Rtsne(obj$pca.rot[,pcs.use], do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else if(dist.use == 'cosine'){
            d = cosine_dist(t(obj$pca.rot[,pcs.use]))
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else {
            d = dist(obj$pca.rot[,pcs.use], method=dist.use)
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        }
    }
    rownames(q) = colnames(obj$data)
    colnames(q) = c('tSNE_1', 'tSNE_2')
    obj$tsne.rot = as.data.frame(q)
    }

    # Cluster cells and run DE tests
    if(length(k) > 0){

        # Save backup singlecell object
	if(do.backup){saveRDS(obj, file=paste0(name, '.obj.rds'))}

	msg( 'Clustering cells', verbose)
	k = k[k < ncol(obj$data)]
	u = paste('Cluster.Infomap.', k, sep='')
	v = run_graph_cluster(data=obj$pca.rot[,pcs.use], k=k)
	obj$meta.data[,u] = v

	if(marker.test != ''){
    	    msg( 'Differential expression', verbose)
            covariates = subset(obj$meta.data, select=c(nGene, Cell_Cycle))
	    obj = set.ident(obj, ident.use=v[,1])
	    print(table(obj$ident))
	    markers = lapply(k, function(ki){
	        obj = set.ident(obj, ident.use=obj$meta.data[,paste0('Cluster.Infomap.', ki)])
	        markers = p.find_all_markers(obj, test.use=marker.test)
	    })
	    names(markers) = k
	}
    }

    if(write_out){

	# Plot TSNE
	png(paste0(name, '.tsne.png'), width=800, height=650)
	plot_tsne(obj, pt.size=1)
	dev.off()

	# Plot clusters
	if(length(k) > 0){
	    pdf(paste0(name, '.clusters.pdf'), width=9, height=9)
	    plot_clusters(obj)
	    dev.off()
	}

	# Marker genes
	if(marker.test != ''){
	    for(ki in names(markers)){write.table(markers[[ki]], file=paste0(name, '.k', ki, '.', marker.test, '.txt'), sep='\t', quote=F)}
	}

	# Save singlecell object
	saveRDS(obj, file=paste0(name, '.obj.rds'))
    }

    return(obj)
}


safeRDS = function(object, file){
    temp = paste0(file, '.temp')
    saveRDS(object, file=temp)
    system(paste('mv', temp, file))
}


merge_obj = function(obj1, obj2, rm_old=FALSE, id1=NULL, id2=NULL){

    # Fix cell names
    if(!is.null(id1)){
	colnames(obj1$counts) = paste(id1, colnames(obj1$counts), sep='.')
	colnames(obj1$data) = paste(id1, colnames(obj1$data), sep='.')
	rownames(obj1$meta.data) = paste(id1, rownames(obj1$meta.data), sep='.')
	rownames(obj1$tsne.rot) = paste(id1, rownames(obj1$tsne.rot), sep='.')
	names(obj1$ident) = paste(id1, names(obj1$ident), sep='.')
    }
    if(!is.null(id2)){
	colnames(obj2$counts) = paste(id2, colnames(obj2$counts), sep='.')
	colnames(obj2$data) = paste(id2, colnames(obj2$data), sep='.')
	rownames(obj2$meta.data) = paste(id2, rownames(obj2$meta.data), sep='.')
	rownames(obj2$tsne.rot) = paste(id2, rownames(obj2$tsne.rot), sep='.')
	names(obj2$ident) = paste(id1, names(obj2$ident), sep='.')
    }

    # Calculate idents
    ident = setNames(c(as.character(obj1$ident), as.character(obj2$ident)), c(colnames(obj1$data), colnames(obj2$data)))
    ident = factor(ident, levels=c(levels(obj1$ident), levels(obj2$ident)))

    # Get metadata
    rows = c(rownames(obj1$meta.data), rownames(obj2$meta.data))
    cols = sort(unique(c(colnames(obj1$meta.data), colnames(obj2$meta.data))))
    meta = matrix(NA, nrow=length(rows), ncol=length(cols))
    rownames(meta) = rows
    colnames(meta) = cols
    meta[rownames(obj1$meta.data), colnames(obj1$meta.data)] = as.matrix(obj1$meta.data)
    meta[rownames(obj2$meta.data), colnames(obj2$meta.data)] = as.matrix(obj2$meta.data)
    meta = as.data.frame(meta)

    # Get tsne coordinates
    tsne = rbind(obj1$tsne.rot[,1:2], obj2$tsne.rot[,1:2])

    # Make counts matrix
    print('Merge counts')
    counts = mem_cbind(list(obj1$counts, obj2$counts))

    # Remove old singlecell objects (save memory)
    if(rm_old == TRUE){rm(obj1); rm(obj2)}

    # Make singlecell object
    obj = make_obj(counts=counts, ming=0, minc=0)
    print('Calculate TPM')
    obj$data = calc_tpm(counts=obj$counts)
    print('Log transform')
    obj$data@x = log2(obj$data@x + 1)
    obj$ident = ident[colnames(obj$data)]

    # Add metadata
    obj$meta.data = as.data.frame(meta)
    obj$tsne.rot = as.data.frame(tsne)

    # Return merged object
    return(obj)
}


merge_objs_hm = function(objs, target='human', ident_fxn=NULL, ortholog_filter='hm'){

    # Merge singlecell objects, using only 1:1 orthologs in human and mouse
    # ---------------------------------------------------------------------
    # Input arguments:
    # - objs = list of singlecell objects ('human' or 'mouse')
    # - target = organism to map genes to ('human' or 'mouse')
    # - ident_fxn = identity function
    # - ortholog_filter:
    #   - 'all' = require orthologs to be found in all singlecell objects
    #   - 'hm'  = require orthologs to be found in human and mouse
    #   - 'none' = no ortholog filter

    # Fix input arguments
    objs = as.list(objs)
    orgs = sapply(objs, function(a) predict_organism(rownames(a$data)))
    print(names(objs))
    print(orgs)

    # Calculate metadata
    cells.use = unname(unlist(sapply(objs, function(a) colnames(a$data))))
    ident.use = unname(unlist(sapply(names(objs), function(a) paste(a, as.character(objs[[a]]$ident), sep='.'))))
    levels(ident.use) = unname(unlist(sapply(names(objs), function(a) paste(a, as.character(levels(objs[[a]]$ident)), sep='.'))))
    org.use = unlist(sapply(1:length(objs), function(i) rep(orgs[[i]], ncol(objs[[i]]$data))))
    dset.use = unlist(sapply(names(objs), function(a) rep(a, ncol(objs[[a]]$data))))
    names(ident.use) = names(org.use) = names(dset.use) = cells.use

    # Define human and mouse gene sets
    h_genes = readLines('/broad/smillie-data/code/csmillie/db/hg19_genes.txt')
    m_genes = readLines('/broad/smillie-data/code/csmillie/db/mm10_genes.txt')

    # Map genes
    ortho_fn = '/broad/smillie-data/code/csmillie/db/hm_orthologs.rds'
    if(file.exists(ortho_fn)){
        print('Reading 1:1 orthologs')
        res = readRDS(ortho_fn)
	h2m = res$h2m
	m2h = res$m2h
    } else {
        print('Calculating 1:1 orthologs')
        h2m = sapply(h_genes, function(gene) {
            mi = map_gene(gene, source='human', target='mouse', do.unlist=F)[[1]]
            hi = map_gene(gene, source='mouse', target='human', do.unlist=F)[[1]]
            if(length(mi) == 1 & length(hi) == 1){
                if(gene == hi){
                    mi
                }
            }
        })
        h2m = sapply(h2m[lengths(h2m) == 1], '[[', 1)
        m2h = setNames(names(h2m), h2m)
	saveRDS(list(h2m=h2m, m2h=m2h), file=ortho_fn)
    }

    # Merge datasets
    print('Merging counts')
    counts = sapply(1:length(objs), function(i) {

        # Get genes to use
	if(orgs[[i]] == 'human'){
	    genes.use = intersect(names(h2m), rownames(objs[[i]]$data))
	} else {
	    genes.use = intersect(names(m2h), rownames(objs[[i]]$data))
	}

	# Subset data
	ci = objs[[i]]$counts[genes.use,]

	# Fix rownames
	if(orgs[[i]] == 'human' & target == 'mouse'){
	    rownames(ci) = unname(h2m[rownames(ci)])
	}
	if(orgs[[i]] == 'mouse' & target == 'human'){
	    rownames(ci) = unname(m2h[rownames(ci)])
	}

	ci

    }, simplify=F)
    counts = sparse_cbind(as.list(counts))
    print(dim(counts))

    # Filter genes
    if(ortholog_filter != 'none'){
        print('Filtering genes')
        if(ortholog_filter == 'all'){
	    groups = dset.use
	} else {
	    groups = org.use
	}
	genes.use = sapply(unique(groups), function(a){
	    j = (groups == a)
	    apply(counts[,j], 1, any)
	})
	counts = counts[apply(genes.use, 1, all),]
	print(dim(counts))
    }

    # Make singlecell object
    obj = make_obj(counts=counts, minc=0, maxc=1e9, ming=0, maxg=1e9, ident_fxn=ident_fxn, verbose=T)
    obj$ident = ident.use[colnames(obj$data)]
    obj$meta.data$organism = org.use[colnames(obj$data)]
    obj$meta.data$dataset = dset.use[colnames(obj$data)]

    obj
}


make_mini = function(obj, num_genes=100, num_cells=100, ident.k=NULL){

    # Select genes and cells
    genes.use = rownames(obj$data)[order(rowMeans(obj$data), decreasing=T)[1:num_genes]]
    cells.use = sample(colnames(obj$data), num_cells)

    # Construct mini singlecell object
    mini = make_obj(obj=obj, genes.use=genes.use, cells.use=cells.use, ming=0, minc=0)

    # Set random identities
    if(!is.null(ident.k)){
        ident = sample(1:ident.k, ncol(mini$data), replace=T)
	mini = set.ident(mini, ident.use=ident)
    }
    return(mini)
}


map_ident = function(obj, old_ident){
    old_ident = as.data.frame(old_ident)
    new_ident = data.frame(ident=rep(NA, ncol(obj$data)), row.names=colnames(obj$data))
    i = intersect(rownames(old_ident), rownames(new_ident))
    new_ident[i,1] = as.character(old_ident[i,1])
    new_ident = structure(as.factor(new_ident[,1]), names=rownames(new_ident))
    return(new_ident)
}

fast_ident = function(obj, ident_map, partial=F){
    ident_map = data.frame(stack(ident_map), row.names=1)
    ident_map = ident_map[levels(obj$ident),,drop=F]
    u = as.character(ident_map[,1])
    ident_map[,1] = ave(u, u, FUN=function(x){if(length(x) > 1){paste(x, 1:length(x), sep='_')} else {x}})
    ident = obj$ident
    if(partial == FALSE){
        levels(ident) = ident_map[as.character(levels(ident)),1]
    } else {
        i = as.character(levels(ident)) %in% rownames(ident_map)
	levels(ident)[i] = ident_map[as.character(levels(ident))[i], 1]
    }
    return(ident)
}

update_signatures = function(obj){

    # Predict host (human or mouse)
    genes = rownames(obj$data)

    # Calculate signatures
    obj$meta.data = cbind(obj$meta.data, score_cells(obj, files='cell_cycle'))
    obj$meta.data = cbind(obj$meta.data, score_cells(obj, files='early'))
    obj$meta.data$nGene = colSums(obj$data > 0)

    return(obj)
}

hclust_ident = function(obj=NULL, data.use=NULL, ident.use=NULL, genes.use=NULL, agg='after', dmethod='euclidean', hmethod='complete'){

    # Input: data.use (features x cells) and ident.use (1 x cells), or singlecell object
    # Output: hclust object

    # Setup input data
    if(is.null(data.use)){
        if(is.null(genes.use)){
	    genes.use = obj$var.genes
	}
	genes.use = intersect(rownames(obj$data), genes.use)
	data.use = obj$data[genes.use,]
    }
    if(is.null(ident.use)){
        ident.use = obj$ident
    }

    # hclust on distance matrix
    if(agg == 'before'){print('Aggregating data')
        print(dim(data.use))
	print(length(ident.use))
        data.use = t(data.frame(aggregate(t(as.matrix(data.use)), list(ident.use), mean), row.names=1))
    }
    print('Calculating distance matrix')
    d = dist(as.data.frame(t(as.matrix(data.use))), method=dmethod)
    if(agg == 'after'){print('Aggregating distances')
        d = data.frame(aggregate(as.matrix(d), list(ident.use), mean), row.names=1)
        d = data.frame(aggregate(t(d), list(ident.use), mean), row.names=1)
    }
    print('Running hclust')
    hclust(as.dist(d), method=hmethod)
}


read_tome_vector = function(tome, name) {
    as.vector(unlist(rhdf5::h5read(tome, name)))
}


read_h5ad_dgCMatrix = function(h5ad, target = "/raw.X") {

    library(Matrix)
    library(rhdf5)

    root <- rhdf5::H5Fopen(h5ad)

    i_path <- paste0(target,"/indices")
    p_path <- paste0(target,"/indptr")
    x_path <- paste0(target,"/data")

    print("Reading indices")
    i <- read_tome_vector(root, i_path)
    print("Reading pointers")
    p <- read_tome_vector(root, p_path)
    print("Reading values")
    x <- read_tome_vector(root, x_path)
    print("Reading observations")
    o <- as.vector(rhdf5::h5read(root, "/obs/_index"))
    print("Reading variables")
    v <- as.vector(rhdf5::h5read(root, "/var/_index"))

    print("Reading dimensions")
    dims <- c(length(v), length(o))

    H5Fclose(root)

    print("Assembling dgCMatrix")
    m <- Matrix::sparseMatrix(i = i,
                              p = p,
                              x = x,
			      dims = dims,
                              index1 = FALSE)

    rownames(m) <- v
    colnames(m) <- o

    return(m)
}


read_h5ad = function(fn, obj=TRUE, counts_regex='layer.*count'){

    library(Matrix)
    library(rhdf5)

    # read metadata
    print('reading metadata')
    meta = h5read(fn, '/obs')
    cats = h5read(fn, '/obs/__categories')
    meta = sapply(names(meta), function(a){
        tryCatch({
        if(a %in% names(cats)){
	    b = factor(as.integer(meta[[a]]))
	    levels(b) = cats[[a]]
	} else {
	    b = meta[[a]]
	}
	as.character(b)
	}, error=function(e) {NULL})
    })
    meta = meta[names(meta) != '__categories']
    meta = as.data.frame(do.call(cbind, meta))
    rownames(meta) = h5read(fn, '/obs/_index')
    print(dim(meta))

    # read coordinates
    print('reading tsne coords from "/obsm/X_umap"')
    coords = tryCatch({t(h5read(fn, '/obsm/X_umap'))}, error=function(e){NULL})

    # read counts
    group = grep(counts_regex, sort(unique(h5ls(fn)$group)), value=T)
    print('reading counts data from h5 index:')
    print(group)
    if(length(group) == 1){
        counts = read_h5ad_dgCMatrix(fn, target=group)
    } else {
        stop('error: found multiple "count" layers')
    }

    if(obj == FALSE){
        return(list(counts=counts, meta=meta, coords=coords))
    } else {

        # make singlecell object
        obj = make_obj(counts=counts, ming=0, minc=0, maxg=1e10, ident_fxn=function(a){'all'})

	# add metadata
	obj$meta.data = cbind(obj$meta.data, meta)

	# add tsne coordinates
	if(!is.null(coords)){
	    obj$tsne.rot = as.data.frame(coords)
	    rownames(obj$tsne.rot) = colnames(obj$data)
	    colnames(obj$tsne.rot) = c('tSNE_1', 'tSNE_2')
	}
	return(obj)
    }
}


