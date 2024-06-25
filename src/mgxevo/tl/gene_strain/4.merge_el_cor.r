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

fl = sprintf('%s/elenta.cor.i%s.tsv',args$outdir,1:50)
df = sapply(fl, fread, simplify=F)
df = do.call(rbind, df)
fwrite(df, 
           sprintf('%s/elenta.cor.tsv', args$outdir), 
           row.names = FALSE)