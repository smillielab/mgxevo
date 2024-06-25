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

parser = ArgumentParser()
parser$add_argument('--outdir')
args = parser$parse_args()

x = fread(sprintf('%s/elenta.cor.tsv', args$outdir))
REF = readLines(sprintf('%s/el.core.txt', args$outdir))
setkey(x, gene)

mu  = mean(x[n.obs>=100][REF]$rho, na.rm=T)
std = sd(x[n.obs>=100][REF]$rho, na.rm=T)
tc = mu - 1.5*std

pred= x[rho>tc][n.obs>=100][pval<5e-3][,unique(gene)]

writeLines(pred,sprintf('%s/el.predicted_genes.txt', args$outdir))