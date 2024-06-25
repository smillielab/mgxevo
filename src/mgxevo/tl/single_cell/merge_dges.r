require(Matrix)
require(optparse)


# ----------------------------------------------------------------------------------
# Combine columns of sparse or dense expression matrices with fast rbind/cbind joins
# ----------------------------------------------------------------------------------
#
# Must specify sample mapping file ("map") or ("name", "path", and "pattern")
# If "path" points to a file, reads as dense matrix with ffread
# Otherwise, reads as sparse matrix with read_mtx(prefix = path)
#
# Optionally:
# - "ming" (min genes per cell) and "minc" (min cells per gene) filters
# - "rename" removes all prefixes from cell barcodes
# - "sparse" saves output as a sparse matrix with writeMM
#
# Usage:
# Rscript merge_dges.r --path /path/to/cellranger/hg19 --out out.txt
# Rscript merge_dges.r --map /path/to/sample/map --out out.mtx --sparse


# ---------------
# Input arguments
# ---------------

option_list = list(make_option('--name', help='prefix added to cell names', default=NULL),
                   make_option('--path', help='sparse matrix path', default=NULL),
		   make_option('--pattern', help='regex pattern for cells to keep', default='.*'),
	           make_option('--map', help='input mapping file (1=prefix, 2=path, 3=pattern)', default='Epi.map.txt'),
		   make_option('--minc', help='cells per gene cutoff', type='integer', default=1),
                   make_option('--ming', help='genes per cell cutoff', type='integer', default=1),
		   make_option('--minr', help='reads per cell cutoff', type='integer', default=1),
                   make_option('--out', help='output file'),
		   make_option('--sep', help='dge delimiter', default='\t'),
		   make_option('--rename', help='use base cell names', action='store_true', default=FALSE),
		   make_option('--rename_sep', help='rename separator', default='\\.'),
		   make_option('--sparse', help='save as sparse matrix', action='store_true', default=FALSE),
		   make_option('--sig', help='filter by gene signature (genes)', default=NULL),
		   make_option('--qmin', help='filter by gene signature (minimum quantile)', type='numeric', default=0),
		   make_option('--transpose', help='transpose matrix?', default=FALSE, action='store_true')

                   )
args = parse_args(OptionParser(option_list=option_list))

require(data.table)
require(Matrix)
require(plyr)
require(tidyverse)
source('~/code/util/mtx.r')
source('~/code/single_cell/scores.r')


# ------------
# Map datasets
# ------------

if(!is.null(args$path)){
    map = matrix(c(args$name, args$path, args$pattern), nrow=1, ncol=3)
} else {
    map = read.table(args$map, stringsAsFactors=F)
}
colnames(map) = c('name', 'path', 'pattern')
map$name[is.na(map$name)] = ''

# -------------
# Read matrices
# -------------

# Read DGE (sparse or text)
read_dge = function(path, prefix='', pattern='.*', ming=1, minc=1, minr=1, rename=FALSE, fix_duplicates=TRUE, transpose=FALSE){

    # Read sparse or text matrix
    if(file.exists(path) & !dir.exists(path)){
        counts = tryCatch({
	    ci = ffread(path, row.names=1, header=TRUE, sep=args$sep)
	    as(as.matrix(ci), 'sparseMatrix')
	}, error=function(e){
	    list(counts=NULL, rmap=NULL)
	})
    } else {
        counts = tryCatch({
	    read_mtx(prefix=path, data='matrix.mtx', rows='features.tsv', cols='barcodes.tsv', filter=TRUE, fix_duplicates=FALSE)
	}, error=function(e){
	    list(counts=NULL, rmap=NULL)
	})
    }

    # Check null
    if(is.null(nrow(counts)) | is.null(ncol(counts))){return(list(counts=NULL, rmap=NULL))}
    if((nrow(counts) == 0) | (ncol(counts) == 0)){return(list(counts=NULL, rmap=NULL))}

    # Transpose
    if(transpose == TRUE){
        counts = t(counts)
    }

    # Filter matrix
    print(dim(counts))
    genes.use = rowSums(counts > 0) >= minc
    cells.use = (grepl(pattern, colnames(counts)) & (colSums(counts > 0) >= ming)) & (colSums(counts) >= minr)
    counts = counts[genes.use, cells.use, drop=F]
    print(dim(counts))

    # Check null
    if(is.null(nrow(counts)) | is.null(ncol(counts))){return(list(counts=NULL, rmap=NULL))}
    if((nrow(counts) == 0) | (ncol(counts) == 0)){return(list(counts=NULL, rmap=NULL))}

    # Fix formatting
    old_names = colnames(counts)
    rownames(counts) = gsub('^mm10_', '', rownames(counts))
    rownames(counts) = gsub('^hg19_', '', rownames(counts))
    colnames(counts) = gsub('\\.1$', '', colnames(counts))
    colnames(counts) = gsub('-1$', '', colnames(counts))
    colnames(counts) = gsub('^\\.', '', colnames(counts))

    # Rename cells
    rmap = list()
    if(rename == TRUE){
        new_names = paste(prefix, gsub(paste0('.*', args$rename_sep), '', colnames(counts)), sep='.')
	colnames(counts) = new_names
	rmap[old_names] = new_names
    } else {
        new_names = paste(prefix, colnames(counts), sep='.')
	colnames(counts) = new_names
	rmap[old_names] = new_names
    }
    colnames(counts) = gsub('^\\.', '', colnames(counts))

    # Fix duplicate genes
    if(fix_duplicates == TRUE){
        print('Fixing duplicates')
        j = which(duplicated(rownames(counts)))
	if(length(j) > 0){
	    counts = as.matrix(counts)
    	    i = match(rownames(counts)[j], rownames(counts))
	    counts[i,] = counts[i,] + counts[j,]
    	    counts = counts[setdiff(1:nrow(counts), j),,drop=F]
	    counts = as(counts, 'sparseMatrix')
	}
    }

    # Filter by gene signature
    if(!is.null(args$sig)){
        sig = intersect(strsplit(args$sig, ',')[[1]], rownames(counts))
	print(paste('Filtering by gene signature:', paste(sig, collapse=',')))
        tpm = log2(calc_tpm(counts=counts) + 1)
	scores = colSums(tpm[sig,,drop=F])
	qmin = quantile(scores, args$qmin)
	cells.use = names(which(scores > qmin))
	counts = counts[,cells.use,drop=F]
	print(dim(counts))
    }

    # Print and return
    cat(paste0('\nRead ', path, ' [', nrow(counts), ' x ', ncol(counts), ']\n'))
    if((args$sparse == TRUE) & (!is.null(dim(counts)))){counts = as(counts, 'sparseMatrix')}
    return(list(counts=counts, rmap=rmap))
}

# size tracking
all_rows = c()
all_cols = c()
all_elem = 0

dges = list()
rmap = list()

for(i in 1:nrow(map)){

    # get dge info
    name = map[i,1]
    path = map[i,2]
    pattern = map[i,3]

    # read and filter dge
    res = read_dge(path, prefix=name, pattern=pattern, ming=args$ming, minc=1, minr=args$minr, rename=args$rename, transpose=args$transpose) # minc=1 for small dges
    dge = res$counts
    print(paste(rownames(dge)[1:5], collapse=', '))
    print(paste(colnames(dge)[1:5], collapse=', '))

    if(is.null(dge) | is.null(ncol(dge))){next}
    rmap[names(res$rmap)] = res$rmap[names(res$rmap)]
    if(!is.null(dim(dge))){dges = c(dges, dge)}

    # size tracking
    all_rows = sort(unique(c(all_rows, rownames(dge))))
    all_cols = sort(unique(c(all_cols, colnames(dge))))
    all_elem = all_elem + length(dge@x)

    print(paste('rows:', length(all_rows), 'cols:', length(all_cols), 'elem:', formatC(all_elem, format='d', big.mark=',')))
}
print(length(dges))


# --------------
# Merge matrices
# --------------

cat('\nMerging DGEs\n')

# sparse_cbind code
# -----------------

# initialize matrix
rows = unique(sort(unlist(sapply(dges, rownames))))
x = Matrix(0, nrow=length(rows), sparse=T)
rownames(x) = rows

# merge with cbind
for(i in 1:length(dges)){
    m = dges[[i]]

    # align rows
    if(length(rows) - nrow(m) > 0){
        n = Matrix(0, nrow=length(rows)-nrow(m), ncol=ncol(m), sparse=T)
        rownames(n) = setdiff(rows, rownames(m))
        m = rbind(m, n)
    }
    m = m[rows,,drop=F]

    # combine data
    x = cbind(x, m)

    # memory management (added this step)
    dges[[i]] = m = 1
}

# remove first column
x = x[,2:ncol(x),drop=F]

# Free up memory
rm(dges)

# Check for replicate cell barcodes
if(max(table(colnames(x))) > 1){
    cat('\nFixing replicates\n')
    cat(paste('\nMerging', ncol(x) - length(unique(colnames(x))), 'replicate cell barcodes'))
    x = as(t(rowsum(t(as.matrix(x)), colnames(x))), 'sparseMatrix')
}


# ------------
# Write output
# ------------

# Filter data
print('Filtering data')
genes.use = rowSums(x > 0) >= args$minc
cells.use = colSums(x > 0) >= args$ming
x = x[genes.use, cells.use]
cat(paste0('\nWriting DGE [', nrow(x), ' x ', ncol(x), ']\n'))

# Write to file
if(args$sparse == TRUE){
    write_mtx(x, prefix=args$out, data='matrix.mtx', rows='features.tsv', cols='barcodes.tsv')
} else {
    cat('\nWriting dense matrix\n')
    require(tibble)
    x = as.data.frame(as.matrix(x)) %>% rownames_to_column('GENE')
    write.table(x, file=args$out, quote=F, sep='\t', row.names=F)
}

# Write cell names
write.table(as.matrix(rmap), file=paste0(args$out, '.rmap.txt'), quote=F, sep='\t', col.names=F)
