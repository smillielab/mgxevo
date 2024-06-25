suppressWarnings({
    suppressMessages({
    library(argparse)

    script_path <- this.path::this.path()
    script_dir <- dirname(script_path)
    mtx <- file.path(script_dir, "/../../ut/mtx.r")
    source(mtx)
    library(data.table)
    })
})

parser = ArgumentParser()
parser$add_argument("-i", type="integer", help="Index")
parser$add_argument("-n", type="integer", help="Number")
parser$add_argument("-g", type="character", help="Gene file")
parser$add_argument("--count_dir", type="character", help="Count directory")
parser$add_argument("--out_dir", type="character", help="Output directory", default='.')

args = parser$parse_args()

i = args$i
n = args$n
g = args$g
count_dir = args$count_dir
out_dir = args$out_dir
out = paste0(out_dir, '/diamond_counts.i', i, '.n', n, '.rds')

chunk = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

genes = sort(unique(readLines(g)))
genes = setdiff(genes, 'Gene')

# map each sample to files
files = list.files(path = count_dir, pattern='\\.counts.txt')
samps = gsub('\\.counts.*', '', gsub('_R[12].*', '', files))
sa2fn = tapply(files, samps, c)
sa2fn = sa2fn[order(names(sa2fn))]

# select samples to use
if(n > 1){
    sa2fn = chunk(sa2fn, n)[[i]]
}

# load sparse matrices
x = sapply(names(sa2fn), function(sa){print(sa)
    xi = sapply(sa2fn[[sa]], function(fn){print(fn)
        u = fread(paste0(count_dir, "/", fn), sep='\t', header=T)
        as.matrix(u[, colnames(u) %in% genes, with=F])
    }, simplify=F)
    xi = colSums(sparse_rbind(xi))
    ci = names(xi)
    xi = matrix(xi, nrow=1)
    rownames(xi) = sa
    colnames(xi) = ci
    xi
}, simplify=F)
x = mem_rbind(x, cleanup=T)

saveRDS(x, file=out)
