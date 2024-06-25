
require(ggplot2)
require(scales)
require(grid)


# Run fast PCA on a singlecell object
run_rpca = function(obj, data=NULL, k, genes.use=NULL, cells.use=NULL, rescale=FALSE, robust=FALSE, scale.max=Inf){

    # This function calculates a PCA on the DGE in obj$scale.data
    # If data, this function will still update the singlecell object
    # If cells.use, you may want to rescale data with rescale=TRUE

    require(rsvd)

    # Get data and subset
    if(is.null(data)){data = obj$scale.data}
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
    data = data[genes.use, cells.use]

    # Calculate PCA
    if(robust == FALSE){
        print(paste('Calculating rpca on [cells x genes] matrix with center =', rescale, 'and scale =', rescale))
	data = t(data)
	if(rescale == TRUE){
	    data = scale(data)
	    data[data < -1*scale.max] = -1*scale.max
	    data[data > scale.max] = scale.max
	}
	pca.obj = rpca(data, center=FALSE, scale=FALSE, retx=TRUE, k=k)
    } else {
        print(dim(data))
        #pca.obj = rrpca(t(data), center=rescale, scale=rescale, retx=TRUE, k=k)
	pca.obj = rrpca(t(data), k=k)
    }
    obj$pca.obj = list(pca.obj)
    obj$pca.rot = data.frame(pca.obj$x)
    obj$pca.x = data.frame(pca.obj$rotation)

    # Return singlecell object
    return(obj)
}

project_pca = function(pca_obj, data, genes.use=NULL, cells.use=NULL, rescale=FALSE, scale_obj=NULL){

    # Project data along the PCs in pca_object

    # Fix PCA object if necessary
    if(length(pca_obj) == 1){pca_obj = pca_obj[[1]]}

    # Subset input data
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
    data = data[genes.use, cells.use]

    # Re-scale input data, if necessary
    if(rescale == TRUE){
        if(is.null(scale_obj)){
	    print('Re-scaling data with center=TRUE, scale=TRUE')
	    data = t(scale(t(data)))
	} else {
	    print('Re-scaling data with parameters in scale_obj')
	    data = t(scale(t(data), center=attr(scale, 'attr:center'), scale=attr(scale, 'attr:scale')))
	}
    } else {print('Assuming data has already been scaled')}

    # Calculate pca.x
    pca.data = t(data) %*% pca_obj$rotation
    return(data.frame(pca.x))
}

# Get significant PCs
get_num_pcs = function(obj, genes.use='', method='karthik', n.cores=1){

    # Use variable genes by default
    if(genes.use == ''){genes.use = obj$var.genes}

    if(method == 'karthik'){
        data = obj$data[genes.use,]
        num_pcs = sig.pcs.perm(t(data), randomized=T, n.cores=n.cores)$r
    }

    if(method == 'jackstraw'){
        num_pcs = jackStraw.permutation.test(obj)$r
    }

    return(num_pcs)
}


# Karthik's significant PC test:
sig.pcs.perm <- function (dat, B = 100, threshold = 0.05,
                                        randomized=F,
                                        verbose=TRUE, seed = NULL,
                                        max.pc=100, n.cores=1,
                                        center=T, scale=T) {

    ptm <- proc.time()
    if(B %% n.cores != 0){stop("Permutations must be an integer multiple of n.cores")}
    cat(sprintf("Scaling input matrix [center=%s, scale=%s]\n", center, scale))
    dat = as.matrix(scale(dat, center=center, scale=scale))
    #dat = as.matrix(t(scale(t(dat), center=center, scale=scale)))
    if (!is.null(seed)) set.seed(seed)
    n <- min(max.pc, ncol(dat))
    m <- nrow(dat)
    print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
    cat(sprintf("Running initial PCA\n"))
    if(randomized){
        require(rsvd)
        uu <- rsvd::rsvd(as.matrix(dat), k=max.pc)
    }else{
        uu <- corpcor::fast.svd(dat, tol = 0)
    }

    ndf <- n - 1
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- matrix(0, nrow = B, ncol = ndf)
    if(verbose==TRUE) message("Estimating number of significant principal components. Permutation: ")

    #permutations
    if(n.cores==1){
        for (i in 1:B) {
            if(verbose==TRUE) cat(paste(i," "))
            dat0 <- t(apply(dat, 1, sample, replace = FALSE))
            if(randomized){
                require(rsvd)
                uu0 <- rsvd::rsvd(as.matrix(dat0), k=max.pc)
            }else{
                uu0 <- corpcor::fast.svd(dat0, tol = 0)
            }
            dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
        }
    }else{
        require(parallel)
        require(foreach)
        require(doParallel)
        cl<-makePSOCKcluster(n.cores, outfile="")
        registerDoParallel(cl, n.cores)
        chunksize = B/n.cores
        vals = split(1:B, ceiling(seq_along(1:B)/chunksize))
        dstat0 = foreach(run.id=1:n.cores, .packages="corpcor", .combine=cbind) %dopar% {
            v = vals[[run.id]]
            #cat(sprintf("Core %s will run perms: %s \n", run.id, paste(v, collapse=",")))
            do.call(rbind, lapply(v, function(i) {
                if(verbose==TRUE) cat(paste(i," "))
                dat0 <- t(apply(dat, 1, sample, replace = FALSE))

                if(randomized){
                    require(rsvd)
                    uu0 <- rsvd(as.matrix(dat0), k=max.pc)
                }else{
                    uu0 <- corpcor::fast.svd(dat0, tol = 0)
                }
                uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
            }))

        }
        cat("\nUnregistering parallel backend..")
        stopCluster(cl)
        registerDoSEQ()
        cat(" done\n");
    }
    p <- rep(1, n)
    for (i in 1:ndf) {
      p[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:ndf) {
      p[i] <- max(p[(i - 1)], p[i])
    }
    r <- sum(p <= threshold)
    y = proc.time() - ptm
    cat(sprintf("\n\n PC permutation test completed. \n %s PCS significant (p<%s, %s bootstraps)\n Runtime: %s s\n ", r,  threshold, B,signif(y[["elapsed"]], 3)))

    return(list(r = r, p = p))
}

# extract and visualise the variance explained by each PC
variance.explained.plot <- function(pca, show_n_pcs=6)
{
  	pdf("variance_explained.pdf", height=7, width=7)
  	show_n_pcs = 10
  	x = head(data.frame(100 * (pca$sdev^2 / sum(pca$sdev^2))), n=show_n_pcs)
  	colnames(x) = "variance"
  	x$pc_index = as.numeric(rownames(x))
  	info("Drawing variance explained plot..")

  	if(max(x) < 10)
  	{
  		max.y = 10
  		increment = 1
  	}else{
  		max.y = max(x)
  		increment = 5
  	}

  	var_plot = ggplot(x, aes(pc_index, variance)) +
  	  	# geom_line(linetype="dashed",  # Dashed line
  	   #            size = 0.5) +       # Thicker line
  	   #  geom_point(shape = 0,         # Hollow squares
  	   #             size = 6)  +
  		geom_bar(stat="identity") +
  	    theme_bw() + xlab("PC") + ylab("% Variance explained") +
  		scale_x_discrete(limits=seq(1,show_n_pcs))
  		#+
  		#scale_y_continuous(breaks=seq(0,max.y,by=increment))
  	print(var_plot)
  	dev.off()
    return (x$variance)
}

# rebecca's significant PC test
sig.pcs <- function(data, binarise.cutoff=0, center=T, scale=F, do.fast=F)
{
    require(BiRewire)
    pc.data.binary=data
    pc.data.binary[pc.data.binary > binarise.cutoff]=1
    if(binarise.cutoff>0){pc.data.binary[pc.data.binary <= binarise.cutoff]=0}
    pc.data.binary.permuted=t(birewire.rewire.bipartite(pc.data.binary,exact=T,verbose = F))
    data.shuffled=apply(data,1,function(x) sample(x,length(x),replace=F))
    pc.data.binary.permuted[pc.data.binary.permuted == 1]=as.numeric(data.shuffled[data.shuffled > 0])
    # scale permuted data set
    pc.data.perm=data.frame(t(scale(t(pc.data.binary.permuted),center=center,scale=scale)))
    # pca on permuted data set

    if (do.fast) {
        pca.perm.obj=fast.prcomp(pc.data.perm,center=F,scale=F)
    } else {
        pca.perm.obj=prcomp(pc.data.perm,center=F,scale=F)
    }
    return(pca.perm.obj)
}

downsample.pca <- function(data, pcs.use=0, num.genes=20, balanced=FALSE, p.threshold=0.001)
{
    scaled = as.matrix( t(scale(t(data), center=TRUE, scale=TRUE)))
    #print(dim(scaled))
    if(pcs.use==0)
    {
        info("Calculating number of PCs significant below p=%s", p.threshold)
        use.pcs = permutationPA(scaled)$r
    }
    pca = prcomp(scaled)
    pc_scores = pca$x
    genes.keep = list()
    for(i in 1:pcs.use) {
        code=paste("PC",i,sep="")
        sx=pc_scores[order(pc_scores[,code]),]
        loaded.genes=rownames(sx)[c(1:num.genes,(nrow(sx)-num.genes):nrow(sx))]
        genes.keep = unlist(c(genes.keep, loaded.genes))
        #print(genes.keep)
    }
    return (data[genes.keep, ])
}

get.loaded.genes <- function(pca, components=1:6, n_genes=20)
{
    rval = NULL
    loadings <- pca$rotation
    for(pc in components){
        #info(sprintf("Getting %s loaded genes for PC-%s ", n_genes, pc))
        l = loadings[, pc]
        top_loadings = as.character(head(names(sort(l, decreasing=TRUE)), n=n_genes))
        bottom_loadings = as.character(tail(names(sort(l, decreasing=TRUE)), n=n_genes))

        v = data.frame(c(top_loadings, bottom_loadings))

        colnames(v) = c(paste("PC-",pc, sep=""))
        #print(v)
        rn = c(1:n_genes, (length(l)-(n_genes-1)):length(l))
        #print(rn)
        rownames(v) = rn
        if(is.null(rval))
        {
            rval=v
        }else
        {
            rval = cbind.data.frame(rval, v)
        }
    }
    return(rval)
}

visualise.pca <- function(pca, groups,
              draw_biplot=T,
							label_samples=T,
							n_components = 4,
							show.loadings=25,
							top_genes = 20)
{


	v = variance.explained.plot(pca, show_n_pcs=8)


	# Show the genes that make up the PCs
	initial_wd = getwd()
	info("PC loadings:")
  loadings <- pca$rotation

  print(loadings[1:5, 1:5])
	loading_dir = "PC_Loaded_Genes"
	dir.create(loading_dir, showWarnings = FALSE)
    setwd(loading_dir)
	top_genes_file = "top_PC_loaded_genes.pdf"
	pdf(top_genes_file, height=11, width=8.5)

	for(pc in 1:n_components){
  		cat(sprintf("Writing top loaded genes for PC-%i ..\n", pc))
  		l = loadings[, pc]
  		top_loadings = data.frame(head(sort(l, decreasing=TRUE), n=top_genes))
  		colnames(top_loadings) = c(paste("PC-",pc, sep=""))

  		grid.table(top_loadings)
  		grid.newpage()

	}
	dev.off()
	bottom_genes_file = "bottom_PC_loaded_genes.pdf"
	pdf(bottom_genes_file, height=11, width=8.5)
	loadings <- pca$rotation
	for(pc in 1:n_components){
		cat(sprintf("Writing bottom loaded genes for PC-%i ..\n", pc))
		l = loadings[, pc]
		top_loadings = data.frame(head(sort(l, decreasing=FALSE), n=top_genes))
		colnames(top_loadings) = c(paste("PC-",pc, sep=""))

		grid.table(top_loadings)
		grid.newpage()

	}
	dev.off()


	# Construct Bi-Plots (eg PC1 vs PC2)
	setwd(initial_wd)
	pc_scatter_dir = "PC_BiPlots"
	dir.create(pc_scatter_dir, showWarnings = FALSE)
	setwd(pc_scatter_dir)
	scores <- data.frame(groups, pca$x[,1:n_components])
	x = c(1:n_components)
	pairs = combn(x, 2)
	print(pairs)
	n_pairs = ncol(pairs)
	cat(sprintf("Generating %i plots.. \n", n_pairs))
	scores$sample_name <- rownames(scores)
	print(head(scores))
	for(p in 1:n_pairs){
		pair = pairs[, p]
		i = pair[1]
		j = pair[2]
		cat(sprintf("Plotting PCs %i vs %i  [Plot %i of %i] ..\n", i, j, p, n_pairs))
		pdf(paste("PC", i, "_PC", j, ".pdf", sep=""), width = 12 , height = 9)

		# regular pc vs pc scatter:
		cat(" Generating regular scatter plot \n")
    n = length(unique(groups))
    p1 = ggplot(data = scores, aes_string(x=paste("PC", i,sep=""), y=paste("PC", j,sep=""), label="sample_name")) +
  			geom_point(aes(fill = factor(groups), shape=21, size=3)) + theme_bw() +
        xlab(sprintf("PC %s   [%s%% variance explained]", i, round(v[i], 1))) +
        ylab(sprintf("PC %s   [%s%% variance explained]", j, round(v[j], 1))) +
  			scale_fill_manual(values=colorRampPalette(brewer.pal(n, "Set1"))(n)) + scale_colour_manual(values=c("black"), guide=FALSE) +
  			scale_shape_identity() +
  	    guides(fill = guide_legend(override.aes = list(shape = 21, size=3), title="Group"), colour=NULL) +
        ggtitle(sprintf("PCA biplot [PC%s vs PC%s]", i, j))
		if(label_samples){
			 p1 = p1 + geom_text(data = scores, vjust=1.5, size=2)
		}
		# ggbiplot fancier scatter:
		if(draw_biplot)
    {
        cat(" Generating biplot\n")
        g <- ggbiplot(pca,
    				obs.scale = 1,
    				var.scale = 1,
    				choices = c(i,j),
    				pc.biplot = TRUE,
    				groups = groups,
    				show.loadings=show.loadings,
    				ellipse = TRUE)
    		g <- g + scale_color_discrete(name = '') + scale_colour_manual(values=colorRampPalette(brewer.pal(n, "Set1"))(n)) +
        scale_shape_identity() +
        guides(fill = guide_legend(override.aes = list(shape = 21, size=3), title="Cell Type"))
    		g <- g + theme(legend.direction = 'horizontal',
    		               legend.position = 'top', legend.title=element_text("Cell Type")) + theme_bw()
        print(g)
    }
		print(p1)

		dev.off()
	}
	setwd(initial_wd)
}

# assumes the input data matrix has samples as columns and
# measurements (genes) as rows
 convert_to_pca_scores <- function(data_matrix, pcs=20){
 	cat(sprintf("Running PCA .. [pcs=%i].. \n", pcs))
 	pca = prcomp(data_matrix, na.action=na.omit, scale=TRUE, center=TRUE)
 	rval = pca$x[, 1:pcs]
 	return (rval)
 }

# modified this to only show top loaded genes. sep27, 2015
#
#  ggbiplot.r
#
#  Copyright 2011 Vincent Q. Vu.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' Biplot for Principal Components using ggplot2
#'
#' @param pcobj           an object returned by prcomp() or princomp()
#' @param choices         which PCs to plot
#' @param scale           covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale       scale factor to apply to observations
#' @param var.scale       scale factor to apply to variables
#' @param pc.biplot       for compatibility with biplot.princomp()
#' @param groups          optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse         draw a normal data ellipse for each group?
#' @param ellipse.prob    size of the ellipse in Normal probability
#' @param labels          optional vector of labels for the observations
#' @param labels.size     size of the text used for the labels
#' @param alpha           alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle          draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param var.axes        draw arrows for the variables?
#' @param varname.size    size of the text for variable names
#' @param varname.adjust  adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev  whether or not to abbreviate the variable names
#'
#' @return                a ggplot2 plot
#' @export
#' @examples
#'   data(wine)
#'   wine.pca <- prcomp(wine, scale. = TRUE)
#'   print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))
#'
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                      obs.scale = 1 - scale, var.scale = scale,
                      groups = NULL, ellipse = FALSE, ellipse.prob = 0.68,
                      labels = NULL, labels.size = 3, alpha = 1,
                      var.axes = TRUE,
                      circle = FALSE, circle.prob = 0.69,
                      varname.size = 3, varname.adjust = 1.5,
                      show.loadings = 20,
                      varname.abbrev = FALSE, ...)
{

  stopifnot(length(choices) == 2)

  # Recover the SVD
 if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
      nobs.factor <- sqrt(pcobj$N)
      d <- pcobj$svd
      u <- predict(pcobj)$x/nobs.factor
      v <- pcobj$scaling
      d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }

  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])

  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)

  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }

  # Scale the radius of the correlation circle so that it corresponds to
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))

  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }

  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs,
                       sprintf('(%0.1f%% explained var.)',
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))

  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }

  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }

  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }

  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)

  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) +
          xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()

  if(var.axes) {
    # Draw circle
    if(circle)
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'),
                         size = 1/2, alpha = 1/3)
    }

    # Draw directions
    # print(dim(df.v))
    # print(head(df.v))
    df.v["variance"] = sqrt(df.v$xvar^2 + df.v$yvar^2)
    df.v = head(df.v[with(df.v, order(-variance)), ], n=show.loadings)
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')),
                   color = muted('red'))
  }

  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups),
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)
    }
  }

  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))

    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }

  #Label the variable axes
  if(var.axes) {
    g <- g +
    geom_text(data = df.v,
              aes(label = varname, x = xvar, y = yvar,
                  angle = angle, hjust = hjust),
              color = 'darkred', size = varname.size)
  }
  #Change the name of the legend for groups
  if(!is.null(groups)) {
    g <- g + scale_color_brewer(name = deparse(substitute(groups)),
                                palette = 'Dark2')
  }

  # TODO: Add a second set of axes

  return(g)
}
