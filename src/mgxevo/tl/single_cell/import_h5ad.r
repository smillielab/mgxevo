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
    parser$add_argument('--prefix', help='prefix (data=*.matrix.mtx, genes=*.genes.csv, cells=*.barcodes.csv, meta=*.meta.csv)', default='')
    parser$add_argument('--out', help='output file (rds)', required=TRUE)
    parser$add_argument('--data', help='matrix (mtx)')
    parser$add_argument('--genes', help='genes (csv)')
    parser$add_argument('--cells', help='cells (csv)')
    parser$add_argument('--meta', help='metadata (tsv)')
    args = parser$parse_args()

}

# -------------------
# fix input arguments
# -------------------

if(args$prefix != ''){
    args$data = paste0(args$prefix, '.matrix.mtx')
    args$genes = paste0(args$prefix, '.genes.csv')
    args$cells = paste0(args$prefix, '.barcodes.csv')
    args$meta = paste0(args$prefix, '.meta.csv')
}

# -----------------------------
# read input data (from scanpy)
# -----------------------------

data = readMM(args$data)
genes = read.csv(args$genes, header=T)[,1]
cells = read.csv(args$cells, header=T)[,1]
meta = read.csv(args$meta, header=T)

if(nrow(data) != length(genes)){
    data = t(data)
}

rownames(data) = genes
colnames(data) = cells


# ------------------------
# create singlecell object
# ------------------------

obj = make_obj(name=args$prefix, dge=data, minc=0, ming=0)
obj$meta.data = cbind(obj$meta.data, meta)


# ----------------
# add pca and umap
# ----------------

pca = paste0(args$prefix, '.pca.csv')
umap = paste0(args$prefix, '.umap.csv')

if(file.exists(pca)){
    pca = read.table(pca, sep=',', header=F)
    rownames(pca) = colnames(obj$data)
    colnames(pca) = paste0('PC_', 1:ncol(pca))
    obj$pca.rot = pca
}

if(file.exists(umap)){
    umap = read.table(umap, sep=',', header=F)
    rownames(umap) = colnames(obj$data)
    colnames(umap) = paste0('tSNE_', 1:ncol(umap))
    obj$tsne.rot = umap
}


# ---------
# save file
# ---------

saveRDS(obj, file=args$out)
