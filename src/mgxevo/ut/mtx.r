library(Matrix)
library(qlcMatrix)

mem_cbind = function(M, cleanup=FALSE){
    require(Matrix)
    
    # Initialize matrix
    x = M[[1]]
    
    # Fast merge
    for(i in 2:length(M)){
        m = M[[i]]
	
	# Add rows to x
	r = setdiff(rownames(m), rownames(x))
	if(length(r) > 0){
	    n = Matrix(0, nrow=length(r), ncol=ncol(x), sparse=T)
	    rownames(n) = r
	    x = rbind(x, n)
	}
	
	# Add rows to m
	r = setdiff(rownames(x), rownames(m))
	if(length(r) > 0){
	    n = Matrix(0, nrow=length(r), ncol=ncol(m), sparse=T)
	    rownames(n) = r
	    m = rbind(m, n)
	}
	
	# Align rows
	m = m[rownames(x),,drop=F]
	
	# Combine data
	x = cbind(x, m)

	# Cleanup data
	if(cleanup){
	    M[[i]] = NA
	}
    }
    
    x
}
sparse_agg = function(x, g, type = "sum"){
    if (type == "sum") {
        x = fac2sparse(g) %*% x
    }
    else if (type == "mean") {
        x = fac2sparse(g) %*% x
        x = x/c(table(g)[rownames(x)])
    }
    return(x)
}

mem_rbind = function(M, cleanup=FALSE){
    require(Matrix)
    
    # Initialize matrix
    x = M[[1]]
    
    # Fast merge
    for(i in 2:length(M)){
        m = M[[i]]
	
	# Add columns to x
	r = setdiff(colnames(m), colnames(x))
	if(length(r) > 0){
	    n = Matrix(0, ncol=length(r), nrow=nrow(x), sparse=T)
	    colnames(n) = r
	    x = cbind(x, n)
	}
	
	# Add columns to m
	r = setdiff(colnames(x), colnames(m))
	if(length(r) > 0){
	    n = Matrix(0, ncol=length(r), nrow=nrow(m), sparse=T)
	    colnames(n) = r
	    m = cbind(m, n)
	}
	
	# Align columns
	m = m[,colnames(x),drop=F]
	
	# Combine data
	x = rbind(x, m)
	
	# Cleanup data
	if(cleanup){
	    M[[i]] = NA
	}
    }
    
    x
}


sparse_cbind = function(M){
    require(Matrix)
    
    # Initialize matrix
    rows = unique(sort(unlist(sapply(M, rownames))))
    x = Matrix(0, nrow=length(rows), sparse=T)
    rownames(x) = rows
    
    # Fast merge with rbind/cbind
    for(m in M){
        
        # Align rows
	if(length(rows) - nrow(m) > 0){
	    n = Matrix(0, nrow=length(rows)-nrow(m), ncol=ncol(m), sparse=T)
	    rownames(n) = setdiff(rows, rownames(m))
	    m = rbind(m, n)
	}
	m = m[rows,,drop=F]
	
	# Combine data
	x = cbind(x, m)
    }

    # Remove first column
    x = x[,2:ncol(x),drop=F]
    x
}


sparse_rbind = function(M){
    require(Matrix)
    
    # Initialize matrix
    cols = unique(sort(unlist(sapply(M, colnames))))
    x = Matrix(0, ncol=length(cols), sparse=T)
    colnames(x) = cols

    # Fast merge with rbind/cbind
    for(m in M){

        # Align columns
	n = Matrix(0, ncol=length(cols)-ncol(m), nrow=nrow(m), sparse=T)
	colnames(n) = setdiff(cols, colnames(m))
	m = cbind(m, n)
	m = m[,cols,drop=F]

	# Combine data
	x = rbind(x, m)
    }

    # Remove first row
    x = x[2:nrow(x),,drop=F]
    x
}

sparse_matsum = function(M){
    require(Matrix)
    
    # Initialize matrix
    rows = unique(sort(unlist(sapply(M, rownames))))
    cols = unique(sort(unlist(sapply(M, colnames))))
    x = Matrix(0, nrow=length(rows), ncol=length(cols), sparse=T)
    rownames(x) = rows
    colnames(x) = cols

    # Fast merge with rbind/cbind
    for(m in M){        
	i = intersect(rownames(x), rownames(m))
	j = intersect(colnames(x), colnames(m))
	x[i,j,drop=F] = x[i,j,drop=F] + m[i,j,drop=F]
    }
    x
}


fast_rnorm = function(x){
    require(wordspace)
    scaleMargins(x, rows=1/rowSums(x))
}


fast_cnorm = function(x){
    require(wordspace)
    scaleMargins(x, cols=1/colSums(x))
}


mtx_filenames = function(prefix, data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv'){
    
    # Get sparse matrix filenames from prefix and patterns
    
    # Get prefix
    if(dir.exists(prefix)){prefix = paste0(prefix, '/')}
    if(!file.exists(prefix)){prefix = paste0(prefix, '.')}
    
    # Fix names
    prefix = gsub('//', '/', prefix)
    #prefix = gsub('\\.\\.', '.', prefix)
    prefix = gsub('/\\.', '/', prefix)
    
    # Get filenames
    if(!file.exists(data)){data = paste0(prefix, data)}
    if(!file.exists(rows)){rows = paste0(prefix, rows)}
    if(!file.exists(cols)){cols = paste0(prefix, cols)}
    
    # Check gunzip extension
    if(!file.exists(data)){data = paste0(data, '.gz')}
    if(!file.exists(rows)){rows = paste0(rows, '.gz')}
    if(!file.exists(cols)){cols = paste0(cols, '.gz')}

    # Remove gunzip extension
    if(!file.exists(data)){data = gsub('\\.gz', '', data)}
    if(!file.exists(rows)){rows = gsub('\\.gz', '', rows)}
    if(!file.exists(cols)){cols = gsub('\\.gz', '', cols)}

    return(list(data=data, rows=rows, cols=cols))
}


read_mtx = function(prefix, data='matrix.mtx', rows='features.tsv', cols='barcodes.tsv', filter=FALSE, fix_duplicates=FALSE){

    # Get filenames
    fns = mtx_filenames(prefix=prefix, data=data, rows=rows, cols=cols)
    print(fns)
    
    # Fix 10X inconsistencies
    if(!all(sapply(fns, file.exists))){
        fns = mtx_filenames(prefix=prefix, data=data, rows='genes.tsv', cols='barcodes.tsv')
    }
    data = fns$data
    rows = fns$rows
    cols = fns$cols
    
    # Read data
    print(paste('Reading', data))
    data = readMM(data)
    print(dim(data))
    rows = read.table(rows, sep='\t', stringsAsFactors=F, quote='', comment.char='')
    
    # Fix 10X inconsistencies
    if(ncol(rows) > 1){
        rows = rows[,sapply(rows, function(a) ! all(a == 'Gene Expression')),drop=F]
	rows = rows[,ncol(rows)]
    } else {
        rows = rows[,1]
    }
    print(length(rows))
    
    cols = read.table(cols, stringsAsFactors=F)
    if(ncol(cols) > 1){
        cols = cols[,ncol(cols)]
    } else {
        cols = cols[,1]
    }
    print(length(cols))
    
    # Filter matrix
    if(filter == TRUE){
        print('read_mtx: minimal filtering (requiring at least 10 genes per cell)')
	j = colSums(data > 0) >= 10
	data = data[,j]
	cols = cols[j]
    }
    
    # Set names
    if(fix_duplicates == TRUE){
        print('read_mtx: fixing duplicate rownames')
	data = as(rowsum(as.matrix(data), rows), 'sparseMatrix')
    } else {
        rownames(data) = rows
    }
    colnames(data) = cols
    
    return(data)
}


read_mtx2 = function(prefix, filter=FALSE, fix_duplicates=FALSE){
    
    # Get filenames
    fns = list.files(pattern=prefix)
    data = grep('matrix', fns, value=T)
    rows = grep('genes|features', fns, value=T)
    cols = grep('barcodes', fns, value=T)
    
    # Read data
    print(paste('Reading', data))
    data = readMM(data)
    if(grepl('csv', rows)){rsep=','} else {rsep='\t'}
    if(grepl('csv', cols)){csep=','} else {csep='\t'}
    rows = read.table(rows, sep=rsep, stringsAsFactors=F, quote='', comment.char='', header=F)
    cols = read.table(cols, sep=csep, stringsAsFactors=F, quote='', comment.char='', header=F)
    
    # Select column with most unique entries
    if(ncol(rows) > 1){
        rows = rows[,which.max(apply(rows, 2, function(a) length(unique(a))))]
    }
    if(ncol(cols) > 1){
        cols = cols[,which.max(apply(cols, 2, function(a) length(unique(a))))]
    }

    # Switch rows and columns if needed
    if(abs(length(rows) - nrow(data)) > (length(cols) - nrow(data))){
        print('Switching rows and columns')
	data = t(data)
    }
    
    # Smart alignment
    if((length(rows) - 1) == nrow(data)){
        rows = rows[2:length(rows)]
    }
    if((length(cols) - 1) == ncol(data)){
        cols = cols[2:length(cols)]
    }    
    
    # Filter matrix
    if(filter == TRUE){
        print('read_mtx: minimal filtering (requiring at least 10 genes per cell)')
	j = colSums(data > 0) >= 10
	data = data[,j]
	cols = cols[j]
    }
    
    # Set names
    if(fix_duplicates == TRUE){
        print('read_mtx: fixing duplicate rownames')
	data = as(rowsum(as.matrix(data), rows), 'sparseMatrix')
    } else {
        rownames(data) = rows
    }
    colnames(data) = cols
    
    return(data)
}


write_mtx = function(x, prefix='.', data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv', temp=FALSE){

    # Get filenames
    if(temp == TRUE){prefix = tempfile(tmpdir='~/tmp')}
    fns = mtx_filenames(prefix=prefix, data=data, rows=rows, cols=cols)
    data = fns$data
    rows = fns$rows
    cols = fns$cols
    
    # Write data
    writeMM(x, data)
    writeLines(rownames(x), rows)
    writeLines(colnames(x), cols)

    # Return filenames
    list(data=data, rows=rows, cols=cols)
}


colMeans_drop0 = function(x){
  n = diff(x@p)
  n[n == 0] = 1
  setNames(colSums(x) / n, 1:ncol(x))
}


rowMeans_drop0 = function(x){
  i = x@i + 1
  n = tabulate(i)
  n[n == 0] = 1
  setNames(rowSums(x) / n, 1:nrow(x))
}


align_mtx = function(x, rows.use=NULL, cols.use=NULL, by_row=TRUE, by_col=TRUE, missing=NA){

    # align a list of matrices by rows and/or columns
    # -----------------------------------------------

    # get row and column names
    if(is.null(rows.use)){
        rows.use = sort(unique(unlist(lapply(x, rownames))))
    }
    if(is.null(cols.use)){
        cols.use = sort(unique(unlist(lapply(x, colnames))))
    }
    
    # check if all matrices are sparse
    sparse = all(sapply(x, function(a) is(a, 'sparseMatrix')))
    
    # align matrices along rows/columns
    y = lapply(x, function(xi){
        if(! by_row){rows.use = rownames(xi)}
	if(! by_col){cols.use = colnames(xi)}
	if(sparse == TRUE){
	    yi = Matrix(missing, nrow=length(rows.use), ncol=length(cols.use), sparse=TRUE)
	    rownames(yi) = rows.use
	    colnames(yi) = cols.use
	    yi[rownames(xi), colnames(xi)] = as(xi, 'sparseMatrix')
	} else {
	    yi = matrix(missing, nrow=length(rows.use), ncol=length(cols.use))
	    rownames(yi) = rows.use
	    colnames(yi) = cols.use
	    yi[rownames(xi), colnames(xi)] = as.matrix(xi)
	}
	yi
    })
    y
}

## AK 03 JAN 2022
# function to pad matricies
# 14 FEB 2023, common cols/rows
pad_mtx  = function(X, rows, cols){
    missing.cols = setdiff(cols,colnames(X))
    missing.rows = setdiff(rows,rownames(X))

    common.cols = intersect(cols, colnames(X))
    common.rows = intersect(rows, rownames(X))
    
    X = X[common.rows,,drop=F]
    X = X[,common.cols,drop=F]
    
    if(length(missing.cols) > 0){
        # pad cols
        yi = Matrix(0, nrow=nrow(X), ncol=length(missing.cols), sparse=TRUE)
        rownames(yi) = rownames(X)
        colnames(yi) = missing.cols
        X = cbind(X,yi)
    }
    if(length(missing.rows) > 0){
        # pad rows
        xi = Matrix(0, nrow=length(missing.rows), ncol=ncol(X), sparse=TRUE)
        rownames(xi) = missing.rows
        colnames(xi) = colnames(X)
        X = rbind(X,xi)
    }
    # align result
    X = X[rows,]
    X = X[,cols]
    return(X)
}
                        
##--------------------------
## Source: https://github.com/saketkc/blog/blob/main/2022-03-10/SparseSpearmanCorrelation2.ipynb
##--------------------------
sparse.ranks <- function(X) {
    if (class(X)[1] != "dgCMatrix") {
        X <- as(object = X, Class = "dgCMatrix")
    }
    non_zeros_per_col <- diff(x = X@p)
    n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
    offsets <- (n_zeros_per_col - 1) / 2
    x <- X@x
    ## split entries to columns
    col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
    ## calculate sparsified ranks and do shifting
    sparsified_ranks <- unlist(x = lapply(X = seq_along(col_lst), 
                                          FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]))
    ## Create template rank matrix
    X.ranks <- X
    X.ranks@x <- sparsified_ranks
    return(X.ranks)
}
sparse_corr <- function(X, Y = NULL, cov = FALSE, method = "pearson") {
 
    if(!(method %in% c("pearson", "spearman"))){
        stop('Method must be either Spearman or Pearson')
    }
    if(method=='pearson'){
        R = qlcMatrix::corSparse(X,Y)
    } else if(method=='spearman'){
        # Get sparsified ranks
        rankX <- sparse.ranks(X)
        if (is.null(Y)){
            # Calculate pearson correlation on rank matrices
            R=(qlcMatrix::corSparse(X=rankX, cov=cov))
        }
        rankY <- sparse.ranks(Y)
        R = (qlcMatrix::corSparse(X = rankX, Y = rankY, cov = cov))
    }
    rownames(R) = colnames(X)
    if(is.null(Y)){
        colnames(R) = colnames(X)
    }
    colnames(R) = colnames(Y)
    # calculate p-values from F-test
    # source: https://stat.ethz.ch/pipermail/r-help/2001-November/016201.html
    P = 1-pf(R^2*(nrow(X)-2)/(1- R^2), 1, nrow(X)-2)
    return(list(r=R, P=P))
}