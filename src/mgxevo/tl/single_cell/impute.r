source('~/code/util/mtx.r')
source('~/code/sge/rsub.r')


impute_magic = function(data, num_pcs=20, k=30, ka=10, eps=1, rescale=99, sparse=TRUE, do.log=FALSE, out=NULL, ret=TRUE, cleanup=TRUE, qsub=FALSE, m=32, t='2:00:00', ...){
    require(data.table)

    # Run magic imputation on TPM or log2TPM (authors use TPM)
    # see: https://github.com/KrishnaswamyLab/magic
    # data = [genes x cells] matrix
    # k = #nn, t = #transitions, ka = autotune, eps = epsilon, rescale = %rescale
    # if is.null(out), then write to temporary file *that gets deleted*

    # check arguments
    if(is.null(out) & ret == FALSE){stop('must specify out with ret=FALSE')}
    if(!is.null(out)){
        if(dir.exists(out)){stop('invalid outfile (is a directory)')}
	if(file.exists(out)){cat('\n*warning* will overwrite outfile\n')}
    }

    # remove zero genes
    data = data[rowSums(data > 0) >= 1,]
    print(dim(data))

    # keep track of cells and filenames
    cells = colnames(data)
    files = c()

    # write data
    print('Writing data')
    if(sparse == FALSE){

        # write data
        tmp = tempfile(tmpdir='~/tmp', fileext='.txt')
	if(is.null(out)){out = tmp}
	files = c(files, tmp)
	fwrite(as.data.frame(t(as.matrix(data))), file=tmp, sep=',', quote=F)

	# magic command
	cmd = paste('python ~/code/single_cell/run_magic.py --txt', tmp)

    } else {

        # write data
	fns = write_mtx(as(t(data), 'sparseMatrix'), temp=TRUE)
	files = c(files, fns$data, fns$rows, fns$cols)
	if(is.null(out)){out = fns$data}

	# magic command
	cmd = paste('python ~/code/single_cell/run_magic.py --mtx', fns$data, '--genes', fns$cols)
    }
    print('Writing data')
    print(files)

    # magic command
    cmd = paste(cmd, '-p', num_pcs, '-k', k, '--ka', ka, '-e', eps, '-r', rescale, '--out', out)

    # submit command
    print('Running magic')
    print(cmd)
    if(qsub == TRUE){rsub(cmd, H='use .tcltk8.6.4', m=m, t=t, ...)} else {system(cmd)}

    # load results
    if(ret == TRUE){
        print('Loading results')
        data = ffread(out, sep=',')
        data = t(data)
        colnames(data) = cells
    } else {
        data = NULL
    }

    # cleanup files
    print('Cleanup')
    if(cleanup == TRUE){for(fn in fns){system(paste0('rm ', fn))}}

    # return results
    return(data)
}


