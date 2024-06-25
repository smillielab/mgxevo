library(argparse)
singlecell()

if(interactive()){

    # set default arguments
    # ---------------------
    args = list()

} else {

    # read input arguments
    # --------------------
    parser = ArgumentParser()
    parser$add_argument('--seur', help='prefix (data=*.matrix.mtx, genes=*.genes.csv, cells=*.barcodes.csv, meta=*.meta.csv)', default='')
    parser$add_argument('--out', help='output file (rds)', required=TRUE)
    args = parser$parse_args()

}


# -----------------------------
# read input data (from seurat)
# -----------------------------
seur = readRDS(args$seur)
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
    seur = data.frame()
}


# ---------
# save file
# ---------

saveRDS(obj, file=args$out)
