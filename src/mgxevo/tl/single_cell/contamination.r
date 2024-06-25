plot_contamination = function(u1, u2, coefs, residuals, lab.use=NULL, lab.fit=NULL, lab.n=NULL, fit.cutoff=2, lab.name='Group', out=NULL){

    require(ggplot2)
    require(ggrepel)
    require(cowplot)

    # log-transform
    l1 = log2(u1 + .5*min(u1[u1 > 0]))
    l2 = log2(u2 + .5*min(u2[u2 > 0]))

    # contamination
    lab.con = names(which(residuals < fit.cutoff))

    # scatterplot data
    d = data.frame(x=l2, y=l1, lab=ifelse(names(l1) %in% lab.fit, names(l1), ''), Type=rep('Other', length(l1)), stringsAsFactors=F)
    d[lab.con, 'Type'] = 'Contamination'
    d[lab.fit, 'Type'] = 'Fit'
    lab.use = intersect(rownames(d), lab.use)
    d[lab.use, 'Type'] = 'Label'
    d[lab.use, 'lab'] = lab.use

    # select subset to label
    if(!is.null(lab.n)){
        lab.sub = sample(lab.fit, min(lab.n, length(lab.fit)))
	d[(!d$lab %in% lab.use & !d$lab %in% lab.sub), 'lab'] = ''
    }

    # rug plot data
    d.rug = data.frame(x=l2[u1 == 0], y=l2[u1 == 0])

    # make plot
    if(!is.null(out)){alpha=.25} else {alpha=1}
    p = ggplot(d, aes(x=x, y=y)) +
        geom_point(aes(colour=Type)) +
   	geom_text_repel(aes(label=lab), size=3, segment.color='grey') +
	#geom_rug(data=d.rug, aes(x=x), sides='t', col='black', alpha=alpha) +
	xlab(paste0('log2(<TPM>) (Non-', lab.name, ')')) +
    	ylab(paste0('log2(<TPM>) (', lab.name, ')')) +
    	scale_colour_manual(values=c('lightcoral', 'black', 'steelblue3', 'lightgray')) +
    	theme_cowplot()

    # add regression lines
    coefs = as.matrix(coefs)
    for(j in 1:ncol(coefs)){
        x0 = (min(l1, na.rm=T) - coefs[1,j] - fit.cutoff)/coefs[2,j]
	x1 = max(l2, na.rm=T)
	d.line = data.frame(x=c(x0, x1))
        d.line$y = coefs[2,j]*d.line$x + coefs[1,j] + fit.cutoff
	if(j == 1){lty = 'longdash'} else {lty = 'dotted'}
	p = p + geom_line(data=d.line, aes(x=x, y=y), lty=lty)
    }

    # save or display plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=2.25, ncol=2.5)
    } else {
        p
    }
}


detect_contamination = function(tpm, idents, samples, anno=NULL, groups=NULL, global_coefs=NULL, fit.n=50, fit.cutoff=2, do.plot=TRUE, lab.use=NULL, prefix='test', n.cores=1){
    # groups = list of groups to use (e.g. Neurons)

    require(MASS)
    source('~/code/single_cell/parallel.r')
    source('~/code/single_cell/regression.r')

    # initialize variables
    if(is.null(anno)){anno = structure(unique(as.character(idents)), names=unique(as.character(idents)))}
    if(any(! idents %in% anno)){stop('! idents %in% anno')}
    if(is.null(groups)){groups = names(anno)}

    # summarize data
    cat('\n\nDetecting ambient contamination\n\n')

    # iterate over groups
    #res = run_parallel(foreach(group=groups) %dopar% {print(group)
    res = sapply(groups, function(group){
        print(group)
        flush.console()

        # output file
        out = paste(prefix, group, 'fit.pdf', sep='.')

        # subset data
	i = idents %in% anno[[group]]
	j = idents %in% setdiff(idents, anno[[group]])
	#j = idents %in% idents

	# sample frequencies
	f = table(as.factor(samples)[i])
	f = as.matrix(f/sum(f))

	# group mean
	u1 = rowMeans(tpm[,i,drop=F])

	# other mean
	u2 = sapply(unique(samples), function(a){
	    rowSums(tpm[,j & (samples == a),drop=F])/sum(samples == a)
	})

	stopifnot(colnames(u2) == colnames(f))
	u2 = (u2 %*% f)[,1]

	# log-transform
	nice_log2 = function(x){y = log2(x); y[is.infinite(y)] = NA; y}
	l1 = nice_log2(u1)
	l2 = nice_log2(u2)

        # fit boundaries
        lo = quantile(l2[u1 == 0], .9, na.rm=T)
        hi = sort(l2, decreasing=T)[100]
        cat(paste0('\n\tLo Cutoff = ', lo, '\n\tHi Cutoff = ', hi))
        exclude = list(c(-Inf, lo), c(hi, Inf))

	# select points for regression
	lab.fit = names(select_points(l2, l1, n=fit.n, dir='down', nbins=10, loess=T, exclude=exclude))
	cat(paste0('\n\tGenes for regression: ', paste(lab.fit, collapse=', ')))

	# robust linear model
	cat('\n\tFitting rlm')
	fit = rlm(l1[lab.fit] ~ l2[lab.fit])
	coefs = as.matrix(coef(fit))
	print(coefs)

	if(!is.null(global_coefs)){coefs = cbind(global_coefs, coefs)}

	# calculate residuals
	residuals = l1 - (coefs[2,1]*l2 + coefs[1,1])
	lab.con = names(which(residuals < fit.cutoff))
	cat(paste0('\n\tLikely contaminants: ', paste(lab.con, collapse=', ')))
	if(do.plot == TRUE){plot_contamination(u1, u2, coefs, residuals, lab.use=lab.use, lab.fit=lab.fit, fit.cutoff=fit.cutoff, out=out, lab.n=20)}

	# update results
	list(u1=u1, u2=u2, fit=fit, coefs=coefs, residuals=residuals, lab.use=lab.use, lab.fit=lab.fit, lab.con=lab.con)
    #}, n.cores = n.cores)
    }, simplify=F)

    names(res) = groups
    cat(paste('\ndetect_contamination: finished\n', names(res), '\n'))
    return(res)
}


full_detect_contamination = function(tpm, groups, idents, samples, anno=NULL, fit.n=50, fit.cutoff=2, do.plot=TRUE, lab.use=NULL, prefix='test', n.cores=1){

    # Find potential contamination in single cell data
    # fit model to "groups" vector
    # run model on "idents" (with optional annotation map)

    # Fit models to cell groups
    cat('\n\nDetect contamination\n\n')
    cat('\nFitting group models\n')
    res.groups = detect_contamination(tpm, idents=groups, samples=samples, anno=NULL, fit.n=fit.n, fit.cutoff=fit.cutoff, do.plot=do.plot, lab.use=lab.use, prefix=prefix, n.cores=n.cores)

    # Average model across groups
    coefs = sapply(res.groups, function(a) a$coefs)
    print(coefs)
    global_coefs = apply(coefs, 1, median)
    cat('\nModel coefficients\n')
    print(coefs)
    cat('\nCombining coefficients\n')
    print(global_coefs)

    # Run average model on groups
    cat('\nFitting group models\n')
    res.groups = detect_contamination(tpm, groups, samples, anno=NULL, fit.n=fit.n, fit.cutoff=fit.cutoff, do.plot=do.plot, lab.use=lab.use, prefix=prefix, global_coefs=global_coefs, n.cores=n.cores)

    # Run average model on idents (with optional annotation map)
    cat('\nFitting ident models\n')
    res.idents = detect_contamination(tpm, idents, samples, anno=anno, fit.n=fit.n, fit.cutoff=fit.cutoff, do.plot=do.plot, lab.use=lab.use, prefix=prefix, global_coefs=global_coefs, n.cores=n.cores)

    # Return data
    return(list(res.groups=res.groups, res.idents=res.idents))
}


add_contamination = function(markers, res){

    # convert residuals to long format
    ires = sapply(res$res.idents, function(a) a$residuals)
    ires = as.data.table(as.data.frame(ires) %>% rownames_to_column('gene') %>% gather('ident', 'contam.res', -gene))

    # merge with marker list
    merge(markers, ires, by=c('ident', 'gene'))
}


run_scrublet = function(counts){

    out = tempfile(pattern='scrublet.', tmpdir='~/tmp', fileext='')
    print(paste('Writing counts:', out))
    write_mtx(t(counts), prefix=out)

    print('Running Scrublet')
    out = paste0(normalizePath(out), '.matrix.mtx')
    cmd = paste0('python ~/code/single_cell/run_scrublet.py --data ', out, ' --out ', out)
    print(cmd)
    system(cmd)

    print('Reading output')
    res = read.table(out, sep='\t')
    rownames(res) = colnames(counts)
    colnames(res) = c('scrublet_score', 'scrublet_prediction')

    print('Cleaning up')
    system(paste0('rm ', out))

    # Return scrublet results
    return(res)
}
