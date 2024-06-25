# This code subsamples N cells per cell type
# It then does PCA and TSNE on the cell subsets
# Then it projects all cells onto these axes

msg( 'Subsampling cells for PCA', verbose)
ident = obj$meta.data$orig.ident
cells.use = as.character(unlist(tapply(colnames(obj$data), list(ident), function(a){sample(a, min(length(a), maxc))})))
msg( sprintf('Selected %d total cells', length(cells.use)))

msg( 'Subsampled PCA', verbose)
if(num_pcs == 0){num_pcs = sig.pcs.perm(t(obj$scale.data[obj$var.genes,cells.use]), randomized=T, n.cores=ncores)$r}
obj$meta.data$num_pcs = num_pcs
msg( sprintf('Found %d significant PCs', num_pcs), verbose)
obj = run_rpca(obj, k=25, genes.use=obj$var.genes, cells.use=cells.use)

msg( sprintf('Subsampled TSNE', seed), verbose)
tsne.rot = Rtsne(obj$pca.rot[cells.use,1:num_pcs], do.fast=T, max_iter=max_iter, verbose=T, perplexity=perplexity)@Y[,1:2]

msg( 'Projecting cells', verbose)
obj$pca.rot = project_pca(obj$pca.obj, obj$scale.data[var_genes,])
new_cells = setdiff(colnames(obj$data), cells.use)
obj$tsne.rot = data.frame(matrix(NA, nrow=ncol(obj$data), 2), row.names=colnames(obj$data))
obj$tsne.rot[cells.use,] = tsne.rot$Y
obj$tsne.rot[new_cells,] = project_tsne(obj$pca.rot[new_cells,1:num_pcs], obj$pca.rot[cells.use,1:num_pcs], obj$tsne.rot[cells.use,], perplexity=perplexity, n.cores=ncores)
colnames(obj$tsne.rot) = c('tSNE_1', 'tSNE_2')
