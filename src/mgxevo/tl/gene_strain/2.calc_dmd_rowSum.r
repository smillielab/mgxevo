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

parser$add_argument('--metadata')
parser$add_argument('--panel')
parser$add_argument('--abd')
parser$add_argument('--dmd_dir')
parser$add_argument('--outdir')
args = parser$parse_args()

#-------------------
# abundance matrices
#-------------------

# subject metadata
meta = data.table(read.metadata.nice(args$metadata))

# load bacterial metadata 
bact = ffread(args$panel, as.dt=T)

# gene counts
# -----------
# ab.g = gene counts matrix [5653 samples x 1987 genes]
# note: contains redundant genes
ab.g = ffread(args$abd, row.names=1, header=T)
i = rowSums(ab.g) >= 1e3
ab.g = ab.g[i,]

# refine panel to only include taxa that survive filtering
bact = bact[apply(bact[,c('dnaG','rpoB','gyrB')], 1, function(a) all(a %in% colnames(ab.g))),]
setkey(bact, genome)

# species counts
# --------------
# ab.b = species counts matrix [5653 samples x 243 species]
# note: non-redundant
ab.b = do.call(cbind, sapply(bact$genome, function(gi) rowMeans(ab.g[,as.character(bact[gi,c('dnaG','rpoB','gyrB')])]), simplify=F))
i = rowSums(ab.b) >= 1e4
ab.b = ab.b[i,]

# tpm normalization
# -----------------
tpm.b = 1e4*ab.b/rowSums(ab.b)
                             
# row sums for DIAMOND
# --------------------                             
i = rownames(tpm.b)
R = setNames(0*(1:len(i)), i)
for(n in 1:50){
    M = readRDS(sprintf('%s/diamond_counts.i%s.n50.rds',args$dmd.dir, n))
    R = R + na.replace(setNames(rowSums(M, na.rm = T)[i], i),0)
}
saveRDS(R, sprintf('%s/rowSums.dmd.rds', args$outdir))
