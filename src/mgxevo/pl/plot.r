script_path <- this.path::this.path()
script_dir <- dirname(script_path)


require(ggplot2)
require(ggrepel)
require(gtools)
require(cowplot)
require(tidyr)


color = paste0(script_dir,'/../tl/single_cell/colors.r')
map_gene = paste0(script_dir,'/../tl/single_cell/map_gene.r')
scores = paste0(script_dir,'/../tl/single_cell/scores.r')

source(color); 
source(map_gene); 
source(scores);

qtrim = function(x, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, rescale=NULL){

    # Trim by value
    x[x < vmin] = vmin
    x[x > vmax] = vmax

    # Trim by quantile
    u = quantile(x, qmin, na.rm=T)
    v = quantile(x, qmax, na.rm=T)
    x[x < u] = u
    x[x > v] = v

    return(x)
}


rescale_vector = function(x, target=c(0,1), f=identity, abs=FALSE){
    # Rescale vector onto target interval
    a = target[[1]]
    b = target[[2]]
    if(min(x) == max(x)){
        rep(min(target), length(x))
    } else {
        if(abs == TRUE){
            x = (x - min(abs(x)))/(max(abs(x)) - min(abs(x)))*(b - a) + a
        } else {
            x = (x - min(x))/(max(x) - min(x))*(b - a) + a
        }
        f(x)
    }
}


subsample_points = function(coords, k, nbins=25, bin_type='size'){

    # Evenly subsample points across multiple axes
    # - coords = data.frame(x1=axis1, x2=axis2, ...)
    # - k = number of points to subsample
    # - nbins = number of bins per axis
    # - bin_type = 'size' (equal size) or 'freq' (equal frequency)
    # For example, splits an xy plot into boxes and samples points from each box
    # Return: TRUE/FALSE vector for each row of coords

    # Divide points into bins
    if(bin_type == 'size'){
        g = interaction(lapply(coords, cut, breaks=nbins))
    } else if(bin_type == 'freq'){
        g = interaction(lapply(coords, cut2, g=nbins))
    } else {
        print("Error: ! bin_type %in% c('size', 'freq')")
    }

    # Subsample points from each bin
    i = as.integer(simple_downsample(1:nrow(coords), groups=g, total_cells=k))
    1:nrow(coords) %in% i
}


load_signature = function(file=NULL, file.regex=NULL, file.cols=NULL, sep=''){
    if(!file.exists(file)){file = paste0('~/aviv/db/markers/', file, '.txt')}
    sig = read.table(file, stringsAsFactors=F, row.names=1, sep=sep, quote='')
    sig = structure(strsplit(sig[,1], ','), names=rownames(sig))
    if(!is.null(file.regex)){file.cols = grep(file.regex, names(sig), value=T)}
    if(!is.null(file.cols)){sig = sig[file.cols]}
    return(sig)
}


resize_jupyter = function(w=10, h=8){
    library(repr)
    options(repr.plot.width=w, repr.plot.height=h)
}

plot_images = function(obj, ..., pt.size=1, titles.use=NULL, do.image=TRUE){
    sink('/dev/null')
    if(is.null(titles.use)){
        titles.use = setNames(names(obj$image), names(obj$image))
    }
    ps = sapply(names(obj$image), function(a){
        if(do.image == FALSE){
	    image.use = NULL
	} else {
	    image.use = obj$image[[a]]
	}
	plot_tsne(obj, image=image.use, coords=obj$impos[[a]], pt.size=pt.size, ...) + ggtitle(titles.use[[a]])
    }, simplify=F)
    legend = get_legend(ps[[1]])
    p = plot_grid(plotlist=lapply(ps, function(a) a + theme(legend.position='none', axis.text.x=element_blank(), axis.text.y=element_blank()) + xlab('') + ylab('')))
    p = plot_grid(p, legend, rel_widths=c(.85, .15))
    sink()
    p
}


plot_tsne = function(obj=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, ident=TRUE, data.use='log2',
                     combine_genes='mean', cells.use=NULL, ymin=0, ymax=1, num_col='auto', pt.size=.75, font.size=11, do.label=T, label.size=5, do.title=TRUE, title.use=NULL, dpal=NULL, cpal=NULL,
	             do.legend=TRUE, legend.title='', share_legend=FALSE, share_legend_rescale=TRUE, legend.size=1, legend_width=.05, legend.cols=NULL, vmin=NA, vmax=NA,
		     na.value='transparent', out=NULL, nrow=1.5, ncol=1.5, return_plotlist=FALSE, ds_cells=NULL, agg.ident=FALSE, do.sort=FALSE, image=NULL, image_ds=10, xzoom=NULL, yzoom=NULL, imborder=50,
		     do.density=FALSE, xlab='Axis 1', ylab='Axis 2', ...){


    # TSNE coordinates
    d = structure(obj$tsne.rot[,1:2], names=c('x', 'y'))
    
    # Cell identities
    if(!is.logical(ident)){
        d$Identity = ident
    } else if(ident & !is.null(obj)){
        d$Identity = obj$ident
    }
        
    # Use alternative coordinates
    if(!is.null(coords)){
        if(is.null(rownames(coords))){
	    rownames(coords) = colnames(obj$data)
	}
	print(head(rownames(d)))
	print(head(rownames(coords)))
	print(length(rownames(d)))
	print(length(rownames(coords)))
	print(length(intersect(rownames(d), rownames(coords))))
	i = intersect(rownames(d), rownames(coords))
	d = d[i,]
	d[rownames(coords),c('x','y')] = coords[,1:2]
    }
        
    # Fix zoom coordinates
    if(!is.null(xzoom)){
        xzoom[1] = max(1, xzoom[1])
	xzoom[2] = min(dim(image)[2], xzoom[2])
    }
    if(!is.null(yzoom)){
        yzoom[1] = max(1, yzoom[1])
	yzoom[2] = min(dim(image)[1], yzoom[2])
    }
        
    # Cell scores
    if(!is.null(cells.use)){cells.use = intersect(cells.use, rownames(d))}
    scores = score_cells(obj=obj, data=data, meta=meta, names=names, regex=regex, files=files, file.cols=file.cols, file.regex=file.regex, top=top, scores=scores,
                         data.use='log2', combine_genes=combine_genes, cells.use=cells.use)

    if(!is.null(scores)){
        i = intersect(rownames(d), rownames(scores))
	ni = c(names(d), names(scores))
	d = cbind.data.frame(d[i,], scores[i,], stringsAsFactors=F)
	colnames(d) = ni
    }
    
    if(agg.ident == TRUE){
        j = setdiff(colnames(d), c('x', 'y', 'Identity'))
        d[,j] = apply(d[,j,drop=F], 2, function(a){
	    ave(a, obj$ident, FUN=function(b){
	        if(is.character(b) | is.factor(b)){names(sort(table(b), dec=T))[[1]]} else{mean(b, na.rm=T)}
	    })
	})
    }

    # Subset cells
    if(is.null(cells.use)){
        cells.use = rownames(d)
    }
    if(!is.null(ds_cells)){
        print(paste('Downsampling', ds_cells, 'cells evenly across xy-grid'))
        i = subsample_points(coords=d[,1:2], k=ds_cells)
	cells.use = cells.use[i]
	print(paste('Selected', length(cells.use), 'cells'))
    }

    if(is.null(cells.use)){cells.use = rownames(d)}
    d = data.frame(d[cells.use,])

    # Initialize plotlist
    cat('\nPlotting:', paste(colnames(subset(d, select=-c(x,y))), collapse=', '), '\n')
    ps = list()

    # Shuffle point order
    d = d[sample(1:nrow(d)),]
    
    # Get limits for shared legend
    if(share_legend == TRUE & share_legend_rescale == TRUE){
        j = names(which(sapply(subset(d, select=-c(x,y)), is.numeric)))
	cat('\nShared limits:', paste(j, collapse=', '), '\n')
        if(is.na(vmin)){vmin = na.omit(min(d[,j]))}
	if(is.na(vmax)){vmax = na.omit(max(d[,j]))}
	cat('> vmin =', vmin, '\n> vmax =', vmax, '\n')
    }

    for(col in setdiff(colnames(d), c('x', 'y'))){

        # plot NAs first
        d = d[c(which(is.na(d[,col])), which(!is.na(d[,col]))),]
	
	# order points by value
	if(do.sort == TRUE){
	    d = d[order(d[,col]),]
	}
	
	if(is.numeric(d[,col])){

	    # Continuous plot

	    # Get colors
	    if(!is.null(cpal)){
	        if(length(cpal) == 1){cont.colors = brewer.pal(9, cpal)} else {cont.colors = cpal}
	    } else {
	        cont.colors = material.heat(50)
	    }
	    d[,col] = qtrim(d[,col], qmin=ymin, qmax=ymax, vmin=vmin, vmax=vmax)

	    p = ggplot(d)
	    
	    if(!is.null(image)){
	        p = p + plot_image(image, xzoom=xzoom, yzoom=yzoom, downsample=image_ds)
	    }
	    
	    p = p +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		scale_colour_gradientn(colours=cont.colors, guide=guide_colourbar(barwidth=.5, title=legend.title), na.value=na.value, limits=c(vmin, vmax)) +
		theme_cowplot(font_size=font.size) +
		xlab(xlab) + ylab(ylab)
	
	} else {

	    # Discrete plot

	    # Get colors
	    if(do.label == TRUE){disc.colors = set2.colors} else {disc.colors = set2.colors}
	    if(!is.null(dpal)){
	        if(is.null(names(dpal))){
		    # use alternative color palette
		    disc.colors = dpal
		} else {
		    # if named palette, then map each name to groups of "similar" colors
		    groups = sapply(names(dpal), function(a) grep(a, levels(d[,col]), value=T))
		    groups = groups[lengths(groups) > 0]
		    dpal = setNames(sapply(1:length(groups), function(i) nice_colors(n=length(groups[[i]]), col=dpal[[i]], type='single')), names(groups))
		    dpal = unlist(lapply(names(dpal), function(a){b = intersect(levels(d[,col]), groups[[a]]); setNames(dpal[[a]], b)}))
		    dpal = dpal[as.character(levels(d[,col]))]
		    disc.colors = dpal
		}
	    }
	    if(! is.factor(d[,col])){d[,col] = factor(d[,col], levels=naturalsort(unique(d[,col])))}

	    p = ggplot(d)

	    if(!is.null(image)){
	        p = p + plot_image(image, xzoom=xzoom, yzoom=yzoom, downsample=image_ds)
	    }
	    
	    p = p +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		theme_cowplot(font_size=font.size) +
		xlab('Axis 1') + ylab('Axis 2') +
		scale_colour_manual(values=disc.colors, na.value=na.value, drop=F) +
		theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
		guides(colour=guide_legend(ncol=legend.cols, title=legend.title))

	    if(do.label == T){
	        t = aggregate(d[,c('x', 'y')], list(d[,col]), median)
		colnames(t) = c('l', 'x', 'y')
		p = p + geom_text_repel(data=t, aes(x=x, y=y, label=l, lineheight=.8), point.padding=NA, size=label.size, family='Helvetica') + theme(legend.position='none')
	    }
	}

	if(do.title == TRUE){
	    if(is.null(title.use)){title = col} else {title = title.use}
	    p = p + ggtitle(title)
	}
	if(do.legend == FALSE){p = p + theme(legend.position='none')}
	if(legend.title == ''){p = p + theme(legend.title=element_blank())}

	if(!is.null(image)){
	    if(!is.null(xzoom)){
	        ixmin = xzoom[1]
		ixmax = xzoom[2]
	    } else {
		ixmin = min(d$x) - imborder
		ixmax = max(d$x) + imborder
	    }
	    if(!is.null(yzoom)){
	        iymin=yzoom[1]
		iymax=yzoom[2]
	    } else {
		iymin = min(d$y) - imborder
		iymax = max(d$y) + imborder
	    }

	    p = p + xlim(c(ixmin, ixmax)) + ylim(c(iymin, iymax))
	}
	
	ps[[col]] = p
    }

    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    
    ps = make_compact(plotlist=ps, num_col=num_col)
    if(length(ps) > 1){
        if(share_legend == TRUE){
            p = share_legend(ps, num_col=num_col, width=legend_width)
        } else {
            p = plot_grid(plotlist=ps, ncol=num_col, align='h')
        }
    }

    if(is.null(out)){
        p
    } else {
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }

    if(return_plotlist == TRUE){return(ps)} else {return(p)}
}


plot_image = function(image, xzoom=NULL, yzoom=NULL, downsample=10){
    if(!is.null(xzoom)){
        xmin = xzoom[1]
	xmax = xzoom[2]
    } else {
        xmin = 1
	xmax = dim(image)[2]
    }
    if(!is.null(yzoom)){
        ymin = yzoom[1]
	ymax = yzoom[2]
    } else {
        ymin = 1
	ymax = dim(image)[1]
    }
    xdown = max(1, as.integer(downsample*(xmax-xmin)/dim(image)[2]))
    ydown = max(1, as.integer(downsample*(ymax-ymax)/dim(image)[1]))
    i = rev(seq(from=ymin, to=ymax, by=xdown))
    j = seq(from=xmin, to=xmax, by=ydown)
    annotation_raster(image[i,j,], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, interpolate=TRUE)
}


make_compact2 = function(plotlist, num_col, labels=TRUE, ticks=TRUE){

    # Make a "plot_grid" plotlist compact by removing axes from interior plots

    # x-axis
    if(length(plotlist) > num_col){
        i = setdiff(1:length(plotlist), rev(1:length(plotlist))[1:min(num_col, length(plotlist))])
	if(labels == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.title.x=element_blank()))}
	if(ticks == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.x=element_blank()))}
    }

    # y-axis
    if(num_col > 1){
        i = setdiff(1:length(plotlist), which(1:length(plotlist) %% num_col == 1))
	if(labels == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.title.y=element_blank()))}
	if(ticks == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.y=element_blank()))}
    }

    return(plotlist)
}


make_compact = function(plotlist, num_col, xlab='b', ylab='l', title=NULL){

    # Make a "plot_grid" plotlist compact by removing axes from interior points

    # Calculate plot indices
    ntot = length(plotlist)
    ncol = min(ntot, num_col)
    indices = list(
        't' = 1:ncol,
	'r' = which(1:ntot %% ncol == 0),
	'b' = rev(1:ntot)[1:ncol],
	'l' = which(1:ntot %% ncol == 1),
	'a' = 1:ntot
    )
    if(ntot == 1){indices[['l']] = 1}

    # Get directions
    xlab = strsplit(xlab, '')[[1]]
    ylab = strsplit(ylab, '')[[1]]
    if(!is.null(title)){title = strsplit(title, '')[[1]]}

    # Fix xlab
    i = setdiff(1:ntot, unlist(indices[xlab]))
    plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.x=element_blank(), axis.title.x=element_blank()))

    # Fix ylab
    i = setdiff(1:ntot, unlist(indices[ylab]))
    plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.y=element_blank(), axis.title.y=element_blank()))

    # Fix title
    if(!is.null(title)){
        i = setdiff(1:ntot, unlist(indices[title]))
        plotlist[i] = lapply(plotlist[i], function(p) p + theme(plot.title=element_blank()))
    }

    return(plotlist)

}


share_legend = function(plotlist, num_col, rel_widths=NULL, width=.1, align='hv'){

    # align = 'hv' makes plots identical sizes (but adds whitespace)

    # Get first legend in plotlist
    i = min(which(sapply(plotlist, function(p) 'guide-box' %in% ggplotGrob(p)$layout$name)))
    cat(paste('\nUsing shared legend:', names(plotlist)[i], '\n'))
    legend = get_legend(plotlist[[i]])

    # Remove all legends
    plotlist = lapply(plotlist, function(p) p + theme(legend.position='none'))

    # Make combined plot
    if(is.null(rel_widths)){rel_widths = rep(1, length(plotlist))}
    p = plot_grid(plotlist=plotlist, ncol=num_col, align=align, rel_widths=rel_widths)
    p = plot_grid(p, legend, ncol=2, rel_widths=c(1-width, width))

    return(p)
}


plot_heatmap = function(obj=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.regex=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL,
                        cells.use=NULL, do.scale=FALSE, scale_method='max', border='black', Rowv=TRUE, Colv=TRUE, labRow=NULL, labCol=NULL, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, out=NA, ...){

    require(NMF)

    # Cell scores
    if(is.null(group_by)){group_by = obj$ident}
    data = get_scores(obj=obj, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.regex=file.regex, file.cols=file.cols, type=type, cells.use=cells.use, group_by=group_by)

    # Scale data
    if(do.scale == TRUE){
        if(scale_method == 'max'){
	    data = as.data.frame(t(t(data)/apply(data, 2, max)))
	} else {
	    data = scale(data)
	}
    }

    # Adjust scale
    data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)

    # Re-order rows
    i = order(apply(data, 1, which.max))
    data = data[i,,drop=F]

    # Plot heatmap
    if(ncol(data) == 1){Colv=NA}
    aheatmap(data, scale='none', border_color=border, Rowv=Rowv, Colv=Colv, labRow=labRow, labCol=labCol, hclustfun='average', filename=out, gp=gpar(lineheight=.8), ...)
}


heatmap2 = function(x, col='nmf', lmat=NULL, lwid=NULL, lhei=NULL, margins=c(5,5), rowSize=1, colSize=1, show.tree=TRUE, ...){
    library(gplots)

    # Make gplots heatmap.2 look like aheatmap

    # Adjust layout
    if(is.null(lmat)){lmat=matrix(c(0,2,3,1,4,0), nrow=2)}
    if(is.null(lwid)){lwid=c(.05,.9,.05)}
    if(is.null(lhei)){lhei=c(.05,.95)}

    # Label sizes
    cexRow = rowSize*(0.2 + 1/log10(nrow(x)))
    cexCol = colSize*(0.2 + 1/log10(ncol(x)))

    # NMF colors
    if(col == 'nmf'){col = colorRampPalette(nmf.colors)(100)} else {col = colorRampPalette(brewer.pal(9, col))(100)}

    # Hide dendrogram
    if(show.tree == FALSE){dendrogram='none'} else {dendrogram='both'}

    # Plot data
    heatmap.2(x, trace='none', lmat=lmat, lhei=lhei, lwid=lwid, col=col, cexRow=cexRow, cexCol=cexCol, dendrogram=dendrogram, ...)
}

ggheatmap = function(data, Rowv='hclust', Colv='hclust', xlab='', ylab='', xsec=FALSE, ysec=FALSE, xstag=FALSE, xstag_space=.15, ystag=FALSE, ystag_space=.15, title='', legend.title='',
                     pal='nmf', do.legend=TRUE, font_size=7, groups=NULL, discrete=FALSE, dpal=set.colors, trans='identity',
                     out=NULL, nrow=1.25, ncol=1.25, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, symm=FALSE, hclust_met='complete', border='#cccccc', replace_na=NA,
		     labRow=NULL, labCol=NULL, ret.order=FALSE, ret.legend=FALSE, pvals=NULL, max_pval=1, pval_border='black', pval_width=.25, na.value='#cccccc', do.label=F, lab.use=NULL, lab.min=NA, lab.offset=0){

    require(ggplot2)
    require(tidyr)
    require(cowplot)
    require(tibble)

    # Scale values
    if(is.null(rownames(data))){rownames(data) = 1:nrow(data)}
    if(is.null(colnames(data))){colnames(data) = 1:ncol(data)}
    data[is.na(data)] = replace_na
    data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)
    
    # Convert to long format
    x = as.data.frame(data) %>% rownames_to_column('row') %>% gather(col, value, -row)
    x$value = as.numeric(x$value)
    
    # Add labels
    if(do.label){
        if(is.null(lab.use)){
	    x$label = x$value
	} else {
	    l = as.data.frame(lab.use) %>% rownames_to_column('row') %>% gather(col, label, -row)
	    x = merge(x, l, by=c('row', 'col'))
	}
    }
    
    # Facet by group
    if(!is.null(groups)){
        # groups is a named list mapping each x or y value to a group
	ri = sum(x$row %in% names(groups))
	ci = sum(x$col %in% names(groups))
	if(ri > ci){
	    x$group = as.factor(groups[x$row])
	} else {
	    x$group = as.factor(groups[x$col])
	}
    }

    # Merge data and p-values
    if(is.null(pvals)){
        x$pval = border
    } else {
        if(ncol(pvals) == 3){
	    colnames(pvals) = c('row', 'col', 'pval')
	    pvals$row = as.character(pvals$row)
	    pvals$col = as.character(pvals$col)
	} else {
	    pvals = as.data.frame(pvals) %>% rownames_to_column('row') %>% gather(col, pval, -row)
	    pvals$pval = as.numeric(pvals$pval)
	}
	if(length(intersect(x$row, pvals$row)) == 0){
	    colnames(pvals) = c('col', 'row', 'pval')
	}
	x = as.data.frame(merge(as.data.table(x), as.data.table(pvals), by=c('row', 'col'), all.x=TRUE))
	x$pval = ifelse(x$pval <= max_pval, pval_border, NA)
	x$pval[is.na(x$pval)] = border
	x = x[order(x$pval),]
    }

    # Order rows
    if(length(Rowv) > 1){rowv = Rowv; Rowv = 'Rowv'} else {rowv = rev(rownames(data))}
    if(length(Colv) > 1){colv = Colv; Colv = 'Colv'} else {colv = colnames(data)}
    if(nrow(data) <= 2){Rowv = 'none'}
    if(ncol(data) <= 2){Colv = 'none'}
    if(Rowv == 'hclust'){
        rowv = rev(rownames(data)[hclust(dist(data), method=hclust_met)$order])
    }
    if(Colv == 'hclust'){
        colv = colnames(data)[hclust(dist(t(data)), method=hclust_met)$order]
    }
    if(Rowv == 'none'){
        rowv = rev(rownames(data))
    }
    if(Colv == 'none'){
        colv = colnames(data)
    }
    if(Rowv == 'min'){
        rowv = rev(rownames(data)[order(apply(data, 1, which.min))])
    }
    if(Rowv == 'max'){
    	i = match(colnames(data)[apply(data, 1, which.max)], colv)
	rowv = rev(rownames(data)[order(i)])
    }
    if(Colv == 'min'){
        i = match(rownames(data)[apply(data, 2, which.min)], rowv)
	colv = rev(colnames(data)[order(i)])
    }
    if(Colv == 'max'){
	i = match(rownames(data)[apply(data, 2, which.max)], rowv)
	colv = rev(colnames(data)[order(i)])
    }
    Rowv = rowv
    Colv = colv

    # Set order of row and column labels
    x$row = factor(x$row, levels=Rowv)
    x$col = factor(x$col, levels=Colv)

    # Get odd/even indices
    if(length(Rowv) > 1){
        r1 = seq(1, length(Rowv), by=2)
        r2 = seq(2, length(Rowv), by=2)
    } else {
        r1 = r2 = 1
    }
    if(length(Colv) > 1){
        c1 = seq(1, length(Colv), by=2)
        c2 = seq(2, length(Colv), by=2)
    } else {
        c1 = c2 = 1
    }

    # Get plot data
    if(length(pal)==1){if(pal == 'nmf'){pal = rev(colorRampPalette(nmf.colors)(101))[10:101]} else {pal = colorRampPalette(brewer.pal(9, pal))(101)}}

    # Plot significant boxes last
    x = x[rev(order(x$pval != 'black')),]
    
    # Add text labels
    if(! 'label' %in% colnames(x)){
        x$label = x$value + lab.offset
        if(!is.na(lab.min)){x$label[x$label < lab.min] = ''}
    }

    if(discrete == TRUE){
        x$value = as.factor(x$value)
	pal = set.colors
    }
    
    # Plot with geom_tile
    p = ggplot(x) +
        geom_tile(aes(x=as.numeric(col), y=as.numeric(row), fill=value), color=x$pval, size=pval_width) +
	labs(x=xlab, y=ylab, title=title, fill=legend.title) +
	theme_cowplot(font_size=font_size) +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.line=element_blank())
    
    if(!is.null(groups)){
        p = p + facet_grid(group ~ ., scales='free', space='free')
    }
    
    # add labels
    if(do.label == TRUE){
        p = p + geom_text(aes(x=as.numeric(col), y=as.numeric(row), label=label))
    }
    
    # Set scale
    if(discrete == TRUE){
        p = p + scale_fill_manual(values=dpal)
    } else {
    if(is.infinite(vmin)){vmin = min(x$value)}
    if(is.infinite(vmax)){vmax = max(x$value)}
    limits = c(vmin, vmax)
    if(symm == TRUE){
        values = c(min(x$value), 0, max(x$value))
	p = p + scale_fill_gradientn(colours=pal, values=scales::rescale(values), na.value=na.value, trans=trans)
    } else {
        p = p + scale_fill_gradientn(colours=pal, limits=limits, na.value=na.value, trans=trans)
    }
    }
    # Secondary x-axis
    if(xsec == FALSE){
        p = p + scale_x_continuous(breaks=1:length(Colv), labels=Colv, expand=c(0,0))
    } else {
	p = p + scale_x_continuous(breaks=c1, labels=Colv[c1], sec.axis=dup_axis(breaks=c2, labels=Colv[c2]), expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=90, hjust=0, vjust=.5))
    }

    # Secondary y-axis
    if(ysec == FALSE){
        p = p + scale_y_continuous(breaks=1:length(Rowv), labels=Rowv, expand=c(0,0))
    } else {
	p = p + scale_y_continuous(breaks=r1, labels=Rowv[r1], sec.axis=dup_axis(breaks=r2, labels=Rowv[r2]), expand=c(0,0))
    }

    # Add axes with ggrepel (use both sides to fit as many labels as possible)
    if(!is.null(labRow)){
        p = p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	if(!is.logical(labRow)){
	    lab.use = intersect(Rowv, labRow)
	    lab.l = ifelse(Rowv %in% lab.use[seq(from=1, to=length(lab.use), by=2)], Rowv, '')
	    lab.r = ifelse(Rowv %in% lab.use[seq(from=2, to=length(lab.use), by=2)], Rowv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none', plot.margin=margin(l=-10, unit='pt'))
	    axis.l = add_ggrepel_axis(plot=p, lab=lab.l, type='left', font_size=font_size)
	    axis.r = add_ggrepel_axis(plot=p, lab=lab.r, type='right', font_size=font_size)
	    p = plot_grid(axis.l, p, axis.r, nrow=1, rel_widths=c(.15,1,.15), align='h', axis='tb')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
    }

    if(!is.null(labCol)){
        p = p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	if(xsec == TRUE & !is.logical(labCol)){
	    lab.use = intersect(Colv, labCol)
	    lab.u = ifelse(Colv %in% lab.use[seq(from=1, to=length(lab.use), by=2)], Colv, '')
	    lab.d = ifelse(Colv %in% lab.use[seq(from=2, to=length(lab.use), by=2)], Colv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none', plot.margin=margin(t=-10, b=-12.5, unit='pt'))
	    axis.u = add_ggrepel_axis(plot=p, lab=lab.u, type='up', font_size=font_size, do.combine=FALSE)
	    axis.d = add_ggrepel_axis(plot=p, lab=lab.d, type='down', font_size=font_size, do.combine=FALSE)
	    p = plot_grid(axis.u, p, axis.d, nrow=3, rel_heights=c(.15,1,.15), align='v', axis='lr')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
	if(xsec == FALSE & !is.logical(labCol)){
	    lab.use = ifelse(Colv %in% intersect(Colv, labCol), Colv, '')
	    legend = get_legend(p)
	    p = p + theme(legend.position='none')
	    axis.d = add_ggrepel_axis(plot=p, lab=lab.use, type='down', font_size=font_size, do.combine=FALSE)
	    p = plot_grid(p, axis.d, nrow=2, rel_heights=c(1, .15), align='v', axis='lr')
	    p = plot_grid(p, legend, nrow=1, rel_widths=c(.925, .075))
	}
    }
    
    if(xstag == TRUE){
        nudge = rep(c(0, 1), length.out=length(Colv))
	p = add_ggrepel_axis(plot=p, lab=Colv, dir='x', type='down', font_size=font_size, force=0, axis.width=xstag_space, nudge=nudge, ret.legend=ret.legend)
    }

    if(ystag == TRUE){
        nudge = rep(c(0, 1), length.out=length(Rowv))
	p = add_ggrepel_axis(plot=p, lab=Rowv, dir='y', type='left', font_size=font_size, force=0, axis.width=ystag_space, nudge=nudge, ret.legend=ret.legend)
    }

    if(do.legend == FALSE){p = p + theme(legend.position='none')}

    # Save plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }

    if(ret.order == FALSE){p} else {list(p=p, Rowv=Rowv, Colv=Colv)}
}

add_ggrepel_axis = function(plot, lab, dir='both', type='left', lab.pos=NULL, lab.lim=NULL, nudge=0, force=1, axis.width=.10, legend.width=.075, font_size=8, ret.legend=FALSE, do.combine=TRUE){

    library(ggrepel)

    # Fix input arguments
    if(is.null(lab.pos)){lab.pos = 1:length(lab)}
    if(is.null(lab.lim)){lab.lim = c(min(lab.pos)-1, max(lab.pos)+1)}

    # Set label directions
    if(type == 'left'){x=0; y=lab.pos; xlim=c(-1,0); ylim=lab.lim; angle=0; nudge_x=-1*nudge; nudge_y=0}
    if(type == 'up'){x=lab.pos; y=0; xlim=lab.lim; ylim=c(0,1); angle=90; nudge_x=0; nudge_y=nudge}
    if(type == 'right'){x=0; y=lab.pos; xlim=c(0,1); ylim=lab.lim; angle=0; nudge_x=nudge; nudge_y=0}
    if(type == 'down'){x=lab.pos; y=0; xlim=lab.lim; ylim=c(-1,0); angle=90; nudge_x=0; nudge_y=-1*nudge}

    # Get data for ggplot
    d = data.frame(x=x, y=y, lab=lab, dir=dir, nudge_y=nudge_y)

    # Make ggrepel axis
    axis = ggplot(d, aes(x=x, y=y, label=lab)) +
           geom_text_repel(min.segment.length=grid::unit(0,'pt'), color='grey30', size=(font_size-1)/.pt, angle=angle, segment.color='#cccccc', segment.size=.15,
	                   direction=dir, nudge_x=nudge_x, nudge_y=nudge_y, force=force) +
	   scale_x_continuous(limits=xlim, expand=c(0,0), breaks=NULL, labels=NULL, name=NULL) +
	   scale_y_continuous(limits=ylim, expand=c(0,0), breaks=NULL, labels=NULL, name=NULL) +
    	   theme(panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0, 'pt'))

    if(do.combine == FALSE){return(axis)}

    # Get plot legend
    legend = get_legend(plot)
    plot = plot + theme(legend.position='none')

    # Combine plots
    if(type == 'left'){
        plot = plot + scale_y_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(l=-12.5, unit='pt')) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	p = plot_grid(axis, plot, nrow=1, align='h', axis='tb', rel_widths=c(axis.width, 1-axis.width))
	if(ret.legend == FALSE){p = plot_grid(legend, p, nrow=1, rel_widths=c(legend.width, 1-legend.width))}
    }
    if(type == 'up'){
        plot = plot + scale_x_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(t=-10, unit='pt')) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	p = plot_grid(axis, plot, nrow=2, align='v', axis='lr', rel_heights=c(axis.width, 1-axis.width))
	if(ret.legend == FALSE){p = plot_grid(legend, p, nrow=1, rel_widths=c(legend.width, 1-legend.width))}
    }
    if(type == 'right'){
        plot = plot + scale_y_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(r=-10, unit='pt')) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	p = plot_grid(plot, axis, nrow=1, align='h', axis='tb', rel_widths=c(1-axis.width, axis.width))
	if(ret.legend == FALSE){p = plot_grid(p, legend, nrow=1, rel_widths=c(1-legend.width, legend.width))}
    }
    if(type == 'down'){
        plot = plot + scale_x_continuous(limits=lab.lim, expand=c(0,0))
        plot = plot + theme(plot.margin=margin(b=-11.5, t=5, unit='pt')) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	p = plot_grid(plot, axis, nrow=2, align='v', axis='lr', rel_heights=c(1-axis.width, axis.width))
	if(ret.legend == FALSE){p = plot_grid(p, legend, nrow=1, rel_widths=c(1-legend.width, legend.width))}
    }

    if(ret.legend == FALSE){p} else {list(p=p, legend=legend)}
}


plot_dots = function(obj=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, files=NULL, file.regex=NULL, file.cols=NULL, top=NULL, groups=NULL, cells.use=NULL,
	             x=NULL, y=NULL, x2=NULL, y2=NULL, # manually control dot sizes and colors
                     dot_size='alpha', dot_color='mu', font_size=8, reorder=FALSE, order_by='alpha',
	             rescale=FALSE, do.title=FALSE, do.legend=TRUE, xlab='', ylab='', out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, max_size=5, drop_zeros=F){

    require(tibble)

    # Fix input arguments
    if(is.null(groups)){groups = obj$ident}
    if(is.null(cells.use)){cells.use = colnames(obj$data)}

    # Get dot sizes and colors (default = alpha and mu)
    x = score_cells(obj=obj, data=data, meta=meta, scores=scores, names=names, regex=regex, files=files, file.regex=file.regex, file.cols=file.cols, top=top, cells.use=cells.use,
        groups=groups, group_stat=dot_size, make.names=F, drop_zeros=drop_zeros)
    y = score_cells(obj=obj, data=data, meta=meta, scores=scores, names=names, regex=regex, files=files, file.regex=file.regex, file.cols=file.cols, top=top, cells.use=cells.use,
        groups=groups, group_stat=dot_color, make.names=F, drop_zeros=drop_zeros)

    # Save levels for ordering
    cells = rownames(x)
    if(reorder == FALSE){
        feats = rev(colnames(x))
    } else {
        feats = colnames(x)[rev(order(apply(x, 2, which.max), -1*apply(x, 2, max)))]
    }

    #if(coord_flip == TRUE){
    #    cells = rev(cells)
    #	feats = rev(feats)
    #}

    # Re-scale dot size and color
    if(rescale == TRUE){
        x = as.data.frame(scale(x, center=F, scale=sapply(x, max)))
	y = as.data.frame(scale(y, center=F, scale=sapply(y, max)))
    }

    # Convert to long format
    d = x %>% rownames_to_column('Group') %>% gather(Feature, Size, -Group)
    y = y %>% rownames_to_column('Group') %>% gather(Feature, Color, -Group)
    d$Color = y$Color

    # Re-order data
    d$Group = factor(d$Group, levels=cells)
    d$Feature = factor(d$Feature, levels=feats)

    # Dot plot
    p = ggplot(d, aes(x=Group, y=Feature, size=Size, fill=Color)) +
        geom_point(color='black', pch=21, stroke=.25) +
	scale_size_area(dot_size, max_size=max_size) +
	scale_fill_gradientn(dot_color, colours=material.heat(8)) +
	theme_cowplot(font_size=font_size) +
	theme(axis.title=element_blank(), panel.grid.major=element_line(colour='black'), axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + panel_border() +
	background_grid(major='xy')

    if(coord_flip == TRUE){p = p + coord_flip()}

    if(!is.null(out)){save_plot(p, file=out, nrow=nrow, ncol=ncol)}

    return(p)
}


plot_dots2 = function(x1=NULL, y1=NULL, x2=NULL, y2=NULL, de=NULL, big_first=FALSE, replace_na=0, na.value='transparent', symm=FALSE, vmin=NULL, vmax=NULL, qmin=0, qmax=1,
                      smin=NULL, smax=NULL, max_size=5,
                      size_title='Fraction of\nexpressing cells', fill_title='Mean non-zero\nexpression\nlog2(TP10K+1)', pval_title='',
		      Rowv='none', Colv='none', pal=rev(brewer.pal(7, 'RdBu')), out=NULL, nrow=1, ncol=1,
		      coord_flip=FALSE, sig_only=TRUE, xlab='', ylab='', font_size=10, pvals=NULL){

    # input arguments
    # x1 = size matrix or list(row, col, size)
    # x2 = size matrix or list(row, col, size)
    # y1 = color matrix or list(row, col, color)
    # y2 = color matrix or list(row, col, color)
    # xanno = annotation track (named list)
    # yanno = annotation track (named list)
    # pvals = data.frame(row, col)

    # convert to long format
    if(!is.null(de)){
        de = as.data.table(de)
        if(sig_only == TRUE){de = de[padjD < .05]}
	x1 = de[,.(col=ident, row=gene, size=alpha)]
	x2 = de[,.(col=ident, row=gene, size=ref_alpha)]
	y1 = de[,.(col=ident, row=gene, color=coefD)]
	y2 = de[,.(col=ident, row=gene, color=coefD)]
    }
    if(!is.null(x1)){
        if(! all(c('row', 'col', 'size') %in% colnames(x1))){
	    x1 = as.data.frame(x1)
	    cnames = rev(colnames(x1))
	    x1 = x1 %>% rownames_to_column('row')
	    x1$row = factor(x1$row, levels=x1$row)
	    x1 = x1 %>% gather(col, size, -row)
	    x1$col = factor(x1$col, levels=cnames)
	}
        x1 = as.data.table(x1)
    }
    if(!is.null(y1)){
        if(! all(c('row', 'col', 'color') %in% colnames(y1))){
	    y1 = as.data.frame(y1)
	    cnames = rev(colnames(y1))
	    y1 = y1 %>% rownames_to_column('row')
	    y1$row = factor(y1$row, levels=y1$row)
	    y1 = y1 %>% gather(col, color, -row)
	    y1$col = factor(y1$col, levels=cnames)
	}
	y1 = as.data.table(y1)
    }
    if(!is.null(x2)){
        if(! all(c('row', 'col', 'size') %in% colnames(x2))){
	    x2 = as.data.frame(x2)
	    cnames = rev(colnames(x2))
	    x2 = x2 %>% rownames_to_column('row')
	    x2$row = factor(x2$row, levels=x2$row)
	    x2 = x2 %>% gather(col, size, -row)
	    x2$col = factor(x2$col, levels=cnames)
	}
	colnames(x2)[colnames(x2) == 'size'] = 'size2'
	if(is.null(x1)){x1 = x2} else {x1 = merge(x1, x2, by=c('row', 'col'), all=T)}
        x2 = as.data.table(x2)
    }
    if(!is.null(y2)){
        if(! all(c('row', 'col', 'color') %in% colnames(y2))){
	    y2 = as.data.frame(y2)
	    cnames = rev(colnames(y2))
	    y2 = y2 %>% rownames_to_column('row')
	    y2$row = factor(y2$row, levels=y2$row)
	    y2 = y2 %>% gather(col, color, -row)
	    y2$col = factor(y2$col, levels=cnames)
	}
	colnames(y2)[colnames(y2) == 'color'] = 'color2'
	if(is.null(y1)){y1 = y2} else {y1 = merge(y1, y2, by=c('row', 'col'), all=T)}
        y2 = as.data.table(y2)
    }

    # merge data
    d = merge(x1, y1, by=c('row', 'col'), all=T)
    for(j in c('size', 'color', 'size2', 'color2')){if(! j %in% colnames(d)){d[,j] = NA}}
    d = as.data.table(d)

    # replace missing values
    d[,'size'][is.na(d[,'size'])] = replace_na
    d[,'size2'][is.na(d[,'size2'])] = replace_na

    # fix plot order
    if(big_first == TRUE & all(c('size', 'size2') %in% colnames(d))){
        d[, order := size <= size2]
	d[order == TRUE, c('size', 'color', 'size2', 'color2') := .(size2, color2, size, color)]
    }
    else{

    }
    d = as.data.table(d)
    print(d)

    # Get row and column orders
    if(length(Rowv) > 1){rowv = Rowv; Rowv = 'Rowv'} else {rowv = levels(as.factor(d$row))}
    if(length(Colv) > 1){colv = rev(Colv); Colv = 'Colv'} else {colv = rev(levels(as.factor(d$col)))}

    if(Rowv == 'auto'){
        u = d[, .SD[which.max(abs(color)), col], row]
	rowv = rev(u[order(factor(u[,V1], levels=colv)), row])
    }

    if(Colv == 'auto'){
        u = d[, .SD[which.max(abs(color)), row], col]
	colv = u[order(factor(u[,V1], levels=rowv)), col]
    }

    d = d[d$row %in% rowv][d$col %in% colv]
    d$row = factor(d$row, levels=rowv)
    d$col = factor(d$col, levels=rev(colv))
    print(d$col)

    # get palette
    if(length(pal) == 1){pal = brewer.pal(7, pal)}
    ci = c(d$color, d$color2)
    si = c(d$size, d$size2)
    if(is.null(vmin)){vmin = min(ci, na.rm=T)}
    if(is.null(vmax)){vmax = max(ci, na.rm=T)}
    if(is.null(smin)){smin = min(si, na.rm=T)}
    if(is.null(smax)){smax = max(si, na.rm=T)}
    d$color = qtrim(d$color, vmin=vmin, vmax=vmax, qmin=qmin, qmax=qmax)
    d$color2 = qtrim(d$color2, vmin=vmin, vmax=vmax, qmin=qmin, qmax=qmax)
    d$size = qtrim(d$size, vmin=smin, vmax=smax)
    d$size2 = qtrim(d$size2, vmin=smin, vmax=smax)
    border1 = ifelse(d$order, '#000000', '#999999')
    border2 = ifelse(d$order, '#999999', '#000000')
    d$sig = 'P > .05'

    # add p-values
    if(!is.null(pvals)){
	colnames(pvals) = c('row', 'col')
        #pvals$row = make.names(pvals$row)
	#pvals$col = make.names(pvals$col)
	pvals$row = as.character(pvals$row)
	pvals$col = as.character(pvals$col)
        setkeyv(d, c('row', 'col'))
        d[pvals[,.(row, col)], 'sig'] = 'P <= .05'
	d$sig = factor(d$sig, levels=c('P <= .05', 'P > .05'))
	p = ggplot(d) + geom_point(aes(x=col, y=row, size=size, fill=as.numeric(color), color=sig), pch=21, stroke=.5) + scale_color_manual(pval_title, values=c('black', 'grey'))
    } else {
        p = ggplot(d) + geom_point(aes(x=col, y=row, size=size, fill=as.numeric(color)), pch=21, color=border1, stroke=.25)
    }

    p = p + theme_cowplot(font_size=font_size) + xlab(xlab) + ylab(ylab) +
	scale_size_area(size_title, max_size=max_size, limits=c(smin, smax)) +
	theme(panel.grid.major=element_line(colour='black'), axis.text.x=element_text(angle=-45, hjust=0, vjust=.5)) + panel_border() + background_grid(major='xy')

    if(any(d$size2 != 0) & any(d$size1 != d$size2)){p = p + geom_point(aes(x=col, y=row, size=size2, fill=as.numeric(color2)), pch=21, color=border2, stroke=.25)}

    if(symm == TRUE){
        values = c(vmin, 0, vmax)
	p = p + scale_fill_gradientn(fill_title, colours=pal, values=scales::rescale(values), na.value=na.value, limits=c(vmin, vmax))
    } else {
        p = p + scale_fill_gradientn(fill_title, colours=pal, na.value=na.value, limits=c(vmin, vmax))
    }
    if(coord_flip == TRUE){print('flip')
	p = p + scale_x_discrete(limits=levels(d$col)) + scale_y_discrete(limits=levels(d$row)) + coord_flip()
    } else {
        p = p + scale_x_discrete(limits=rev(levels(d$col))) + scale_y_discrete(limits=rev(levels(d$row)))
    }
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }

    p
}


plot_dots3 = function(obj, genes.use, ident.use=NULL, order='auto', Colv=NULL, drop_zeros=T, cells.use=NULL, ...){

    # fix input arguments
    if(is.null(ident.use)){ident.use = obj$ident}

    # expression stats
    alpha = score_cells(obj, genes.use, groups=ident.use, group_stat='alpha', drop_zeros=drop_zeros, cells.use=cells.use)
    mu = score_cells(obj, genes.use, groups=ident.use, group_stat='mu', drop_zeros=drop_zeros, cells.use=cells.use)
    mean = score_cells(obj, genes.use, groups=ident.use, group_stat='mean', drop_zeros=drop_zeros, cells.use=cells.use)

    # order genes
    if(is.null(Colv) & order == 'auto'){
        Colv = order(apply(mean, 2, which.max), -apply(mean, 2, max))
	Colv = colnames(mean)[Colv]
    }
    print(Colv)

    # make dotplot
    plot_dots2(alpha, mu, alpha, mu, big_first=T, Colv=Colv, ...)
}


plot_violin = function(obj=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL,
                       regex=NULL, files=NULL, file.regex=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL, color_by=NULL, pt.size=.25,
	               do.facet=FALSE, facet_by=NULL, facet_genes=NULL, facet_formula=NULL, facet_scales='free_y', qmin=0, qmax=1, vmin=NULL, vmax=NULL, do.scale=FALSE, num_col='auto', resort=NULL,
                       ident=TRUE, cells.use=NULL, do.title=TRUE, title.use=NULL, do.legend=FALSE, xlab='', ylab='log2(TP10K+1)', out=NULL, nrow=1.5, ncol=1.5, legend.title='Group',
		       coord_flip=FALSE, alpha=1, order=NULL, font_size=10, do.pvals=FALSE, combine_genes='scale2', ret.data=FALSE, ds_cells=NULL,
		       mar=unit(c(6,6,6,6),'pt')){

    library(ggbeeswarm)
    # Initialize plots
    ps = list()

    # Get facet formula
    if(is.null(facet_formula)){
    if(do.facet == FALSE){
        facet_formula = ifelse(is.null(facet_by), '~ .', 'Facet ~ .')
    } else {
        do.title = FALSE
        facet_formula = ifelse(is.null(facet_by), 'Feature ~ .', 'Feature ~ Facet + .')
    }
    } else {facet_formula = as.formula(facet_formula)}

    # Fix input arguments
    if(is.null(group_by)){group_by = obj$ident}
    if(is.null(color_by)){color_by = group_by}
    if(is.null(facet_by)){facet_by = rep('', ncol(obj$data))}
    if(is.null(cells.use)){cells.use = colnames(obj$data)}

    # Plot data
    d = data.frame(Group=as.factor(group_by), Color=as.factor(color_by), Facet=as.factor(facet_by), row.names=colnames(obj$data))
    d = d[cells.use,,drop=F]
    if(coord_flip == TRUE){d$Group = factor(d$Group, levels=rev(levels(d$Group)))}

    # Downsample cells
    if(!is.null(ds_cells)){
        print(paste('Downsampling', ds_cells, 'cells'))
	print(table(d$Group))
        cells.use = simple_downsample(cells=rownames(d), groups=d$Group, total_cells=ds_cells)
	d = d[cells.use,,drop=F]
	print(table(d$Group))
    }

    # Cell scores
    scores = score_cells(obj=obj, data=data, meta=meta, scores=scores, names=names, regex=regex, files=files, top=top, file.regex=file.regex, file.cols=file.cols, cells.use=cells.use,
                         combine_genes=combine_genes)
    scores = scores[cells.use,,drop=F]
    names = colnames(scores)
    d = cbind.data.frame(d, scores)

    # Fix NAs
    d = d[!is.na(d$Facet),,drop=F]

    # Scale data
    j = sapply(d, function(a) !is.factor(a))
    d[,j] = apply(d[,j,drop=F], 2, function(a) qtrim(a, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax))

    # Facet data
    if(do.facet == TRUE){
	d = gather(d, Feature, Value, -Group, -Color, -Facet)
	d$Feature = factor(d$Feature, levels=names)
    }

    if(!is.null(facet_genes)){
        d$Facet = facet_genes[d$Feature]
    }

    # Calculate p-values
    d$Group = droplevels(d$Group)
    if(do.pvals == TRUE){
        pvals = do.call(rbind, lapply(levels(d$Feature), function(fi) {do.call(rbind,
        lapply(levels(d$Group), function(gi) {do.call(rbind,
	    lapply(setdiff(levels(d$Facet), levels(d$Facet)[[1]]), function(hi){
                u = d[(d$Feature == fi & d$Group == gi) & d$Facet == levels(d$Facet)[[1]], 'Value']
	        v = d[(d$Feature == fi & d$Group == gi) & d$Facet == hi, 'Value']
                c(fi, gi, hi, max(v), tryCatch({wilcox.test(u,v,use='pairwise.complete.obs')$p.value}, error=function(e){1}))
	    })
	)})
    )}))

    pvals = data.frame(pvals, stringsAsFactors=F)
    colnames(pvals) = c('Feature', 'Group', 'Facet', 'Value', 'Pval')
    pvals$Value = as.numeric(pvals$Value)
    pvals$Pval = p.adjust(as.numeric(pvals$Pval), 'fdr')
    pvals$Pval[is.na(pvals$Pval)] = 1
    pvals$Label = ifelse(pvals$Pval <= 1e-3, '***', ifelse(pvals$Pval <= 1e-2, '**', ifelse(pvals$Pval <= .05, '*', '')))
    pvals$Color = pvals$Group
    pvals$Facet = factor(pvals$Facet, levels=levels(d$Facet))
    }

    # Violin plots
    print(dim(d))
    for(col in setdiff(colnames(d), c('Group', 'Color', 'Facet', 'Value'))){

        # Convert to numeric?
	if(!is.factor(d[,col]) & !is.character(d[,col])){d[,col] = as.numeric(d[,col])}

        # Facet if necessary
        if(do.facet == TRUE){
	    if(!is.null(resort)){
	        i = d$Feature == levels(d$Feature)[[1]]
	        group_order = names(sort(tapply(d$Value[i], d$Group[i], resort)))
		d$Group = factor(d$Group, levels=group_order)
	    }
	    p = ggplot(data=d, aes_string(x='Group', y='Value', fill='Color'))
	} else {
	    if(!is.null(resort)){
	        group_order = names(sort(tapply(d$Group, d[,col], resort)))
		d$Group = factor(d$Group, levels=group_order)
	    }
	    p = ggplot(data=d, aes_string(x='Group', y=col, fill='Color'))
	}

	# Make violin plot
	p = p +
	    #geom_point(position=position_jitterdodge(dodge.width=0.6, jitter.width=1), size=pt.size, show.legend=F) +
	    geom_quasirandom(position=position_dodge(), size=pt.size, show.legend=F, method='pseudorandom') +
	    geom_violin(scale='width', alpha=alpha, size=.25) +
	    scale_fill_manual(values=set.colors) + theme_cowplot(font_size=font_size) +
	    xlab(xlab) + ylab(ylab) + labs(fill=legend.title) +
	    stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom='crossbar', width=.5, show.legend=F, size=.25) +
	    theme(axis.text.x = element_text(angle = -45, hjust = 0))

	if(do.title == TRUE){
	    if(is.null(title.use)){
	        p = p + ggtitle(col)
	    } else {
	        p = p + ggtitle(title.use)
	    }
	}

	if(do.legend == FALSE){
	    p = p + guides(fill=FALSE)
	}

	if(coord_flip == TRUE){
	    p = p + coord_flip()
	}

	if(facet_formula != '~ .'){
	    p = p + facet_grid(as.formula(facet_formula), scales=facet_scales)
	}

	if(do.pvals == TRUE){
	    p = p + geom_text(data=pvals, aes(x=Group, y=Value + .01*max(Value), label=Label))
	}

	ps[[col]] = p
    }
    print(length(ps))

    if(length(ps) > 1){

        if(num_col == 'auto'){
            num_col = ceiling(sqrt(length(ps)))
        }
	if(do.facet == TRUE){
            ps = make_compact(ps, num_col=num_col)
	}
	p = plot_grid(plotlist=ps, ncol=num_col)
        p = p + theme(plot.margin=mar)
    } else {
        p = ps[[1]]
    }

    if(!is.null(out)){
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
    p
}


plot_feats = function(...){plot_tsne(...)}

plot_pcs = function(obj, pcs=NULL, ...){
    if(is.null(pcs)){pcs = 1:obj$meta.data$num_pcs[[1]]}
    pcs = obj$pca.rot[,pcs]
    plot_tsne(obj, scores=pcs, ...)
}

plot_clusters = function(obj, ...){
    clusters = sort(grep('Cluster', colnames(obj$meta.data), value=T))
    clusters = sapply(obj$meta.data[,clusters], function(a){droplevels(as.factor(a))})
    plot_tsne(obj, scores=clusters, ...)
}


matrix_barplot = function(data, group_by=NULL, pvals=NULL, xlab='', ylab='Frequency', value='mean', error='se', legend.title='Groups', colors=NULL, pos='dodge', border=NA,
                          out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, sig_only=F, do.facet=F, facet_formula='group ~ .', do.pie=FALSE){

    # Plot barplot of [M x N] matrix
    # x-axis = matrix columns (e.g. cell types)
    # y-axis = matrix values (e.g. frequencies)
    # fill = matrix rows (e.g. samples) or groups (e.g. conditions)

    # Arguments:
    # group.by = the group of each row
    # pvals = [G x N] matrix of p-values for each group and column
    # error = sd, se, or none

    require(tidyverse)
    require(data.table)
    require(ggsignif)

    # Groups (default = rows)
    if(is.null(group_by)){group_by = rownames(data)}
    if(nlevels(group_by) == 0){group_by = factor(group_by, levels=group_by)}

    # Select significant comparisons
    if(sig_only == TRUE){
        j = apply(pvals, 2, min) <= .05
	if(sum(j) == 0){return(NULL)}
	data = data[,j,drop=F]
	pvals = pvals[,j,drop=F]
    }

    # Construct input data
    names = colnames(data)
    data = data.frame(group=group_by, data)
    group_levels = levels(group_by)
    if(ncol(data) > 2){
        colnames(data)[2:ncol(data)] = names
    }
    data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'group')))
    
    # Value function
    if(value == 'mean'){vf = mean} else if(value == 'median'){vf = median} else {stop()}

    # Error function
    se = function(x, na.rm=T){if(length(x) == 1){0} else {sd(x, na.rm=na.rm)/sqrt(length(x))}}
    if(error == 'sd'){ef = sd} else if(error == 'se'){ef = se} else {ef = function(x, ...){0}}

    # Estimate error bars
    data = data[,.(u=vf(y, na.rm=T), s=ef(y, na.rm=T)),.(group, x)]
    
    # Add p-values 1
    if(!is.null(pvals)){
        pvals = as.data.frame(pvals) %>% rownames_to_column('group') %>% gather(x, pval, -group) %>% as.data.table()
	setkeyv(data, c('x', 'group'))
	setkeyv(pvals, c('x', 'group'))
	data = merge(data, pvals, all=T)
	data$lab1 = ifelse(data$pval <= .001, '***', ifelse(data$pval <= .01, '**', ifelse(data$pval <= .05, '*', '')))
    }

    if(coord_flip == TRUE){names = rev(names); group_levels=rev(group_levels)}
    data$x = factor(data$x, levels=names)
    data$group = factor(data$group, levels=group_levels)

    # Get colors
    if(is.null(colors)){colors = disc.colors(length(group_levels))}

    # Plot data
    if(pos == 'stack'){
        p = ggplot(data[!is.na(group)]) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity', na.rm=T)
	if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', width=.25)}
    } else {
        pos = position_dodge(.9)
        p = ggplot(data[!is.na(group)]) + geom_bar(aes(x=x, y=u, fill=group), colour=border, size=.25, stat='identity', na.rm=T, position=pos)
	if(error %in% c('sd', 'se')){p = p + geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', position=pos, width=.25)}
    }

    p = p +
        scale_fill_manual(values=colors, name=legend.title) + xlab(xlab) + ylab(ylab) +
        scale_color_manual('', values=c('#000000', '#999999', '#cccccc'), guide='none')

    # Facet wrap
    if(do.facet == TRUE){
        p = p + facet_grid(as.formula(facet_formula), scales='free')
    }
    p = p + theme_cowplot()

    dy = max(data$u + data$s, na.rm=T)*.01
    if(coord_flip == FALSE){
        p = p + theme(axis.text.x = element_text(angle = -45, hjust = 0))
	if(!is.null(pvals)){p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group), hjust='center', vjust=0, size=5, angle=0, position=pos)}
    } else {
        p = p + coord_flip()
	if(!is.null(pvals)){p = p + geom_text(aes(x=x, y=u+s+dy, label=lab1, group=group), hjust='center', vjust=1, size=5, angle=90, position=pos)}
    }

    # Save plot
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}


matrix_piechart = function(data, group_by=NULL, pvals=NULL, pcut=0.05, value='mean', do.all=FALSE, font_size=8, pal=NULL, num_row=NULL, num_col=NULL, border=NA, out=NULL, nrow=1.5, ncol=1.5, reorder=FALSE){

    # Plot pieplot of [M x N] matrix
    # x-axis = matrix columns (e.g. cell types)
    # y-axis = matrix values (e.g. frequencies)
    # fill = matrix rows (e.g. samples) or groups (e.g. conditions)
    # group.by = the group of each row

    require(tidyverse)
    require(data.table)

    # Groups (default = rows)
    if(is.null(group_by)){group_by = rownames(data); if(reorder == FALSE){group_by = factor(group_by, levels=rownames(data))}}
    if(nlevels(group_by) == 0){group_by = factor(group_by)}
    
    # Construct input data
    names = colnames(data)
    if(value == 'mean'){vf = mean} else if(value == 'sum'){vf = sum} else {stop('error: invalid "value" argument')}
    data = data.frame(aggregate(data, list(group_by), vf, na.rm=T), row.names=1)
    if(do.all == TRUE){data$All = rowSums(data); names=c(names, 'All')}
    data = data %>% rownames_to_column('group')
    data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'group')))
    data$x = factor(data$x, levels=names)
    data$group = factor(data$group, levels=levels(group_by))
    
    # Get palette
    if(is.null(pal)){pal = set.colors}

    # Add p-values
    if(!is.null(pvals)){
        pvals = as.data.frame(pvals) %>% gather(x, pval) %>% as.data.table()
	setkey(data, x)
	setkey(pvals, x)
        data = merge(data, pvals, all=T)
	data$color = factor(ifelse(data$pval <= pcut, paste0('P < ', pcut), paste0('P > ', pcut)))
    } else {
        data$color = '1'
    }

    # Plot data
    p = ggplot(data[!is.na(group)]) +
        geom_bar(aes(x=1, y=y, fill=group, color=color), stat='identity', position=position_fill()) +
	facet_wrap(. ~ x, nrow=num_row, ncol=num_col) +
	coord_polar('y', start=pi/2, direction=1) +
	scale_y_continuous(expand = c(0,0)) +
	scale_fill_manual('', values=pal) + scale_color_manual('', values=c('#000000', '#cccccc')) +
        theme_cowplot(font_size=font_size) +
	theme(axis.text=element_blank(), axis.ticks=element_blank(), line=element_blank(), strip.background=element_blank()) +
        xlab('') + ylab('') + labs(fill='')

    if(length(unique(data$color)) == 1){p = p + guides(color=FALSE)}

    # Save plot
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}


plot_volcano = function(fcs, pvals, labels=NULL, color_by=NULL, facet_by=NULL, lab.x=c(-1, 1), max_pval=.05, lab.n=0, lab.p=0, lab.size=5, do.repel=TRUE, pval_floor=-Inf,
                        font.size=11, legend.title='Groups', out=NULL, nrow=1.5, ncol=1.5, xlab='Fold change', ylab='-log10(pval)', palette=NULL, ret.labs=FALSE, do.log=T){

    # Volcano plot with automatic labeling:
    # best p-values in region lab.x (n = lab.n)
    # best p-values globally (n = lab.p)

    # Make plot data
    if(is.null(color_by)){color_by = rep('Group', length(fcs))}
    if(is.null(facet_by)){facet_by = rep('Group', length(fcs))}
    pvals[pvals < pval_floor] = pval_floor
    if(do.log == TRUE){pvals = -log10(pvals)}
    data = data.frame(x=fcs, y=pvals, Color=color_by, Facet=facet_by, stringsAsFactors=F)
    data$Facet = as.factor(data$Facet)

    if(!is.null(labels)){
        for(facet in unique(data$Facet)){

	    # Select facet data
	    i = which(data$Facet == facet & 10**(-data$y) <= max_pval)
	    di = data[i,]

	    # Label points < lab.x
	    breaks = seq(from=min(di$x, na.rm=T), to=lab.x[[1]], length.out=10)
	    groups = cut(di$x, breaks=breaks, include.lowest=TRUE)
	    j1 = as.numeric(simple_downsample(cells=1:nrow(di), groups=groups, ngene=di$y, total_cells=as.integer(lab.n/2)))

	    # Label points > lab.x
	    breaks = seq(from=lab.x[[2]], to=max(di$x, na.rm=T), length.out=10)
	    groups = cut(di$x, breaks=breaks, include.lowest=TRUE)
	    j2 = as.numeric(simple_downsample(cells=1:nrow(di), groups=groups, ngene=di$y, total_cells=as.integer(lab.n/2)))

	    # Label best global p-values
	    j3 = which(order(-1*di$y) <= lab.p)

	    # Set labels
	    j = sort(unique(c(j1, j2, j3)))
	    data[i[j], 'labels'] = TRUE
	}
    } else {
        data$labels = FALSE
    }

    # Set labels
    data$labels = ifelse(data$labels, labels, '')

    # Make volcano plot
    p = ggplot(data, aes(x=x, y=y)) +
        geom_point(aes(colour=Color)) +
	theme_cowplot(font_size=font.size) +
	xlab(xlab) +
	ylab(ylab)

    # Add labels
    if(do.repel == TRUE){
        p = p + geom_text_repel(aes(label=labels), size=lab.size, segment.color='#cccccc')
    } else {
        p = p + geom_text(aes(label=labels), size=lab.size)
    }

    # Add color scale
    if(is.numeric(data$Color)){
        if(is.null(palette)){palette = colorRampPalette(rev(brewer.pal(7, 'RdBu')))(100)}
        p = p + scale_colour_distiller(palette=palette)
    } else {
        if(is.null(palette)){palette = desat(set.colors, .5)}
        p = p + scale_colour_manual(values=palette)
    }
    if(all(color_by == 'Group')){p = p + theme(legend.position='none')}

    # Facet wrap
    if(nlevels(data$Facet) > 1){
        print('Faceting')
        p = p + facet_wrap(~ Facet, scales='free') +
	    theme(axis.line = element_line(colour = "black"),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank(),
	    strip.background = element_blank())

    }

    # Save to file
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    if(ret.labs){as.character(na.omit(data$labels))} else {p}
}


plot_volcanos = function(markers, color_by=NULL, outdir='volcano'){

    # Create output directory
    if(!file.exists(outdir)){dir.create(outdir)}
    if(!is.null(color_by) & (nrow(markers) != length(color_by))){stop()}

    # Split markers by contrast
    m = split(markers, markers$contrast)
    q = split(color_by, markers$contrast)

    # Make each volcano plot
    for(name in names(m)){
        out = gsub(':', '.', paste0(outdir, '/', name, '.volcano.png'))
	print(c(name, out))
        plot_volcano(m[[name]]$coefD, m[[name]]$pvalD, labels=m[[name]]$gene, color_by=q[[name]], lab.x=c(-.5, .5), max_pval=.05, lab.n=100, lab.p=20, lab.size=3, out=out, nrow=2, ncol=2)
    }
}


simple_scatter = function(x, y, lab=NA, sig=FALSE, col=NULL, col.title='', size=NULL, size.title='', lab.use=NULL, lab.sig=FALSE, lab.near=0,  lab.n=0, lab.g=0, groups=NULL, lab.size=4, lab.type='up',
                          palette=NULL, legend.fontsize=10, border=F, edges=NULL, na.value='#cccccc', do.label=F,
                          xlab=NULL, ylab=NULL, out=NULL, nrow=1, ncol=1, min_size=1, max_size=3, xskip=c(0,0), yskip=c(0,0), xlim=NULL, ylim=NULL, alpha=1, unlab.grey=FALSE, auto_color='multi',
			  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, xmin.n=1, xmax.n=1, ymin.n=1, ymax.n=1, cmin=-Inf, cmax=Inf,
			  legend.names=c('', 'lab.n', 'lab.use', 'lab.near'), do.sort=FALSE, max.overlaps=10, lab.fontface = "italic"){

    require(ggplot2)
    require(ggrepel)
    require(cowplot)
    require(wordspace)

    x = as.numeric(x)
    y = as.numeric(y)
    if(is.null(col)){col = rep('', length(x))}
    if(is.null(sig)){sig = rep(NA, length(x))}
    if(is.null(size)){size = rep(1, length(x))}
    if(is.character(col)){col = as.factor(col)}
    
    if(!is.null(edges)){colnames(edges) = c('x', 'y', 'xend', 'yend')}
    
    # adjust range
    xmin = max(na.omit(sort(x))[[xmin.n]], xmin)
    xmax = min(na.omit(sort(x, decreasing=T))[[xmax.n]], xmax)
    ymin = max(na.omit(sort(y))[[ymin.n]], ymin)
    ymax = min(na.omit(sort(y, decreasing=T))[[ymax.n]], ymax)
    x[x < xmin] = xmin
    x[x > xmax] = xmax
    y[y < ymin] = ymin
    y[y > ymax] = ymax

    if(!is.factor(col)){
        col = pmax(cmin, col)
        col = pmin(cmax, col)
    }
    
    d = data.frame(x=x, y=y, lab=lab, col=col, size=size, flag='', sig=sig, stringsAsFactors=FALSE)
    i.lab = !is.na(d$lab)
    di = d[i.lab,]

    if(do.sort){d = d[order(d$col),]}

    if(is.null(xlab)){xlab = deparse(substitute(x)); if((length(xlab) > 5) || (substr(xlab, 1, 2) == 'c(')){xlab=''}}
    if(is.null(ylab)){ylab = deparse(substitute(y)); if((length(ylab) > 5) || (substr(ylab, 1, 2) == 'c(')){ylab=''}}

    if(lab.n > 0 | lab.g > 0){

        i = c()

	if('up' %in% lab.type | 'down' %in% lab.type){

            # get breaks
	    if(!is.null(groups)){groups.use = groups} else {
	        breaks = seq(from=min(di$x, na.rm=T), to=max(di$x, na.rm=T), length.out=min(20, lab.n))
	        groups.use = cut(di$x, breaks=breaks, include.lowest=TRUE)
		groups.use[xskip[[1]] < di$x & di$x < xskip[[2]]] = NA
	    }

	    # get cells
	    if('up' %in% lab.type){
	        ngene = ifelse(yskip[[1]] < di$y & di$y < yskip[[2]], NA, di$y)
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=ngene, total_cells=lab.n)))
		i = c(i, order(-1*di$y)[1:lab.g])
	    }
	    if('down' %in% lab.type){
	    	ngene = ifelse(yskip[[1]] < di$y & di$y < yskip[[2]], NA, -1*di$y)
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=ngene, total_cells=lab.n)))
		i = c(i, order(di$y)[1:lab.g])
	    }
	}

	if('right' %in% lab.type | 'left' %in% lab.type){

	    # get breaks
	    if(!is.null(groups)){groups.use = groups} else {
	        breaks = seq(from=min(di$y, na.rm=T), to=max(di$y, na.rm=T), length.out=min(20, lab.n))
	        groups.use = cut(di$y, breaks=breaks, include.lowest=TRUE)
		groups.use[yskip[[1]] < di$y & di$y < yskip[[2]]] = NA
	    }

	    # get cells
	    if('right' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=di$x, total_cells=lab.n)))
		i = c(i, order(-1*di$x)[1:lab.g])
	    }

	    if('left' %in% lab.type){
	        i = c(i, as.numeric(simple_downsample(cells=1:nrow(di), groups=groups.use, ngene=-1*di$x, total_cells=lab.n)))
		i = c(i, order(di$x)[1:lab.g])
	    }
	}

	if('random' %in% lab.type){

	    i = c(i, sample(1:nrow(di), lab.n))
	}

	di[unique(i), 'flag'] = 'lab.n'
    }
    d[i.lab,] = di

    if(!is.null(lab.use)){

        # label neighbors
	if(lab.near > 0){
	    u = as.matrix(d[, c('x','y')])
	    v = as.matrix(d[lab %in% lab.use, c('x','y')])
	    i = unique(sort(unlist(apply(dist.matrix(u,v,skip.missing=TRUE,method='euclidean'), 2, order)[1:(lab.near + 1),])))
	    d$flag[i] = 'lab.near'
	}

        # label points
        d$flag[d$lab %in% lab.use] = 'lab.use'
    }

    if(lab.sig == TRUE){d$flag[d$sig == TRUE] = 'lab.use'; d$sig = FALSE}

    d$lab[d$flag == ''] = ''
    d = d[order(d$flag),]
    d$flag = factor(d$flag, levels=c('', 'lab.n', 'lab.use', 'lab.near'), ordered=T)

    if(auto_color == 'none'){d$col = ''}
    if(auto_color == 'bw' & all(col == '')){d$col = ifelse(d$flag == '', '', 'lab.n')}
    if(auto_color == 'multi' & all(col == '')){d$col = d$flag}

    # plot labeled points last
    d = d[order(d$flag),]
    i = is.na(d$col) | d$col == ''
    d = rbind(d[i,], d[!i,])

    if(unlab.grey == TRUE){d[d$lab == '',]$col = ''}

    # hack: fix legend names
    if(!is.numeric(d$col)){
        d$col = as.factor(d$col)
    	#levels(d$col) = legend.names
    }

    p = ggplot(d)

    if(!is.null(edges)){p = p + geom_segment(data=edges, aes(x=x,y=y,xend=xend,yend=yend), color='black')}
    
    if(border == FALSE){
        if(length(unique(size)) == 1){
	    p = p + geom_point(aes(x=x, y=y, col=col), size=d$size, alpha=alpha)
	} else {
            p = p + geom_point(aes(x=x, y=y, col=col, size=size), alpha=alpha)
	}
    } else {
        if(length(unique(size)) == 1){
	    p = p + geom_point(aes(x=x, y=y, fill=col, stroke=ifelse(sig == TRUE, stroke, 0)), size=size, alpha=alpha, pch=21, colour=border)
	} else {
	    p = p + geom_point(aes(x=x, y=y, fill=col, size=size, stroke=ifelse(sig == TRUE, stroke, 0)), alpha=alpha, pch=21, colour=border)
	}
    }
    p = p + geom_text_repel(data=d[(d$lab != '') | (1:nrow(d) %in% sample(1:nrow(d), min(nrow(d), 1000))),], aes(x=x, y=y, label=lab), size=lab.size,fontface=lab.fontface, segment.color='grey', max.overlaps=max.overlaps) +
	theme_cowplot() + xlab(xlab) + ylab(ylab) + theme(legend.title=element_text(size=legend.fontsize), legend.text=element_text(size=legend.fontsize))

    if(all(size == 1)){p = p + scale_size(guide = 'none', range=c(min_size, max_size))} else {p = p + scale_size(name=size.title, range=c(min_size, max_size))}

    if(!is.null(xlim)){p = p + xlim(xlim)}
    if(!is.null(ylim)){p = p + ylim(ylim)}

    # label colors on plot
    if(do.label){
        t = aggregate(d[,c('x', 'y')], list(d[,'col']), median)
        colnames(t) = c('l', 'x', 'y')
        p = p + geom_text_repel(data=t, aes(x=x, y=y, label=l, lineheight=.8), point.padding=NA, size=lab.size, family='Helvetica', fontface=lab.fontface) + theme(legend.position='none')
    }

    if(border != FALSE){
        if(!is.null(palette)){
            if(!is.numeric(d$col)){
	        p = p + scale_fill_manual(name=col.title, values=palette, na.value=na.value, breaks=levels(as.factor(d$col)))
	    } else {
	        p = p + scale_fill_gradientn(name=col.title, colors=palette, na.value=na.value)
	    }
        } else {
            if(!is.numeric(d$col)){
	        p = p + scale_fill_manual(name=col.title, values=c('lightgrey', 'black', 'red', 'pink')) + theme(legend.position='none')
	    } else {
	        p = p + scale_fill_gradientn(name=col.title, colors=material.heat(100), na.value=na.value)
	    }
        }
    } else {
        if(!is.null(palette)){
	    if(!is.numeric(d$col)){
                p = p + scale_color_manual(name=col.title, values=palette, na.value=na.value, breaks=levels(as.factor(d$col)))
	    } else {
	        p = p + scale_color_gradientn(name=col.title, colors=palette, na.value=na.value)
	    }
        } else {
            if(!is.numeric(d$col)){
	        p = p + scale_color_manual(name=col.title, values=c('lightgrey', 'black', 'red', 'pink')) + theme(legend.position='none')
	    } else {
	        p = p + scale_color_gradientn(name=col.title, colors=material.heat(100), na.value=na.value)
	    }
        }
    }

    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}

    return(p)
}


simple_volcano = function(mi, x='coefD', y='padjD', lab='gene', xskip=c(-1,1), xlab='Fold change\n(discrete coefficient)', ylab='-log10(Adjusted p-value)', ...){
    simple_scatter(mi[[x]], -log10(mi[[y]]), lab=mi[[lab]], xskip=xskip, xlab=xlab, ylab=ylab, yskip=c(0, -log10(.05)), ...) + geom_abline(slope=0, intercept=-log10(.05), linetype='dashed')
}


