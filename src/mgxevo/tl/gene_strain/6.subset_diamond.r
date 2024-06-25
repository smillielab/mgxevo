suppressWarnings({
    suppressMessages({
    library(argparse)

    script_path <- this.path::this.path()
    script_dir <- dirname(script_path)
        
    mtx <- file.path(script_dir, "/../../ut/mtx.r")
    treetool <- paste0(script_dir,'/../tree/treetools.r')
    util <- paste0(script_dir,'/../../ut/util.r')
        
    source(mtx)
    library(data.table)
    source(treetool)
    source(util)
    })
})

parser$add_argument('--meta')
parser$add_argument('--gene')
parser$add_argument('--dmd_dir')
parser$add_argument('--outdir')
parser$add_argument('--outfile')
args = parser$parse_args()


genes = readLines(args$gene)
ofn2  = sprintf('%s/%s', args$outdir, args$outfile)

m = read.metadata.nice(args$metadata)
rows = m$sample
cols = genes

# initialise matrix
M = Matrix(0, nrow=length(rows), ncol=length(cols), sparse=TRUE)
rownames(M) = rows
colnames(M) = cols

# iterate through DIAMOND files
for(n in 1:50){
    M =  M + pad_mtx(readRDS(sprintf('%s/diamond_counts.i%s.n50.rds',args$dmd_dir, n)), rows = rows, cols = cols)
}
saveRDS(M, file = ofn2)



