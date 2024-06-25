

mean_cv_loess = function(data, num_genes=1500, use_bins=TRUE, num_bins=20, window_size=100, do.plot=FALSE, invert=FALSE, use_var=FALSE){

    # calculate mean and cv
    u = apply(data, 1, mean)
    v = apply(data, 1, var)
    i = u > 0 & v > 0
    u = u[i]
    v = v[i]
    if(use_var == TRUE){
        cv = v
    } else {
        cv = sqrt(v)/u
    }

    # fit loess curve
    l = loess(log(cv) ~ log(u), family='symmetric')
    d = log(cv) - l$fitted

    # get variable genes
    if(use_bins == TRUE){

        # select variable genes from equal frequency bins
        library(Hmisc)
	k = as.integer(num_genes/num_bins)
	if(invert == FALSE){
	    var_genes = as.character(unlist(tapply(d, cut2(u, g=num_bins), function(a) names(sort(a, decreasing=T)[1:k]))))
	} else {
	    var_genes = as.character(unlist(tapply(d, cut2(u, g=num_bins), function(a) names(sort(a, decreasing=F)[1:k]))))
	}

    } else {

        # select variable genes with rolling z-scores
        library(zoo)

	# re-order by mean expression
	D = d[order(u)]

	# rolling z-scores (only use unique values)
	ru = rollapply(D, window_size, function(a) mean(unique(a)), partial=T)
	rs = rollapply(D, window_size, function(a) sd(unique(a)), partial=T)
	rz = structure((D - ru)/rs, names=names(D))

	# select variable genes
	var_genes = names(sort(rz, decreasing=T)[1:num_genes])
    }

    # plot results
    if(do.plot == TRUE){
	colors = c('#cccccc', 'black')[as.integer(names(u) %in% var_genes) + 1]
        plot(log(u), log(cv), pch=16, col=colors, xlab='log(mean)', ylab='log(cv)')
	lines(l$x[order(u)], l$fitted[order(u)], col='red', lw=2)
    }

    return(var_genes)
}


get_var_genes = function(data, ident=NULL, use_var=FALSE, genes.use=NULL, method='loess', num_genes=1500, min_cells=5, min_ident=25, do.plot=F, prefix=NULL, do.flatten=T, n.cores=1, ...){

    source('~/code/single_cell/parallel.r')
    
    # Fix ident
    if(is.null(ident)){ident = rep(1, ncol(counts))}
    ident = as.factor(as.character(ident))

    # Filter ident
    u = sort(table(ident))
    i = names(which(u <= min_ident))
    if(sum(u[i]) >= min_ident){
        new_name = 'Merge'
    } else {
        new_name = names(which.min(u[u >= min_ident]))
    }
    levels(ident)[levels(ident) %in% i] = new_name
    print('Calculating var genes across')
    print(table(ident))
    
    # Subset idents
    levels.use = names(sort(table(ident), dec=T))[1:min(50, length(unique(ident)))]
    print(table(ident)[levels.use])

    # Subset genes
    if(is.null(genes.use)){
        genes.use = rownames(data)
    } else {
        genes.use = intersect(rownames(data), genes.use)
	print(paste('Subsetting from', nrow(data), 'to', length(genes.use), 'genes'))
    }
    data = data[genes.use,]

    # var_genes = sapply(levels(ident), function(i){print(i)
    var_genes = sapply(levels.use, function(i){print(i)

            # Start plotting device
	    if(!is.null(prefix)){png(paste(prefix, i, 'png', sep='.'), w=1000, h=800)}

	    # Subsample data
	    data = data[,ident == i]
    	    genes.use = rowSums(data > 0) >= min_cells
	    data = data[genes.use,]
	    print(dim(data))

	    if(method == 'loess'){
	        vi = mean_cv_loess(data, use_var=use_var, num_genes=num_genes, do.plot=do.plot, ...)
	    }
	    if(method == 'adam'){
	        source('~/dev/adam/rna_seq/r/samples.r')
		source('~/dev/adam/rna_seq/r/util.r')
	        vi = get.variable.genes(data, do.plot=do.plot, ...)
		i = (!is.na(var_genes[,4])) & (var_genes[,4] <= .05)
		vi = vi[i,1]
	    }
	    if(method == 'karthik'){
	        vi = meanCVfit(data, diffCV.num_genes=num_genes, do.plot=do.plot, ...)
	    }

	    # Stop plotting device
	    if(!is.null(prefix)){dev.off()}

	    return(vi)
    })

    if(do.flatten == TRUE){
        a = sort(table(as.character(unlist(var_genes))), decreasing=T)
	num_genes = min(length(a), num_genes)
	k = a[num_genes]
	u = names(a)[a > k]
	v = sample(names(a)[a == k], num_genes - length(u))
	var_genes = c(u,v)
    }

    return(var_genes)
}


meanCVfit = function(count.data, reads.use=FALSE, do.text=FALSE, diffCV.cutoff=NULL, diffCV.num_genes=NULL, do.spike=FALSE, main.use=NULL, prefix=NULL, do.plot=FALSE, ret.diffCV=FALSE){

    require(MASS)

    # Empirical mean, var and CV
    mean_emp = apply(count.data, 1, mean)
    var_emp = apply(count.data, 1, var)
    cv_emp = sqrt(var_emp) / mean_emp

    # NB sampling
    a=colSums(count.data)
    size_factor =  a/ mean(a)
    fit=fitdistr(size_factor, "Gamma")
    if (do.spike) spike.genes=grep("^ERCC", rownames(count.data), value=TRUE)
    print(fit)

    if(do.plot==T){
    par(mfrow=c(2,2))
    if (!reads.use){
        hist(size_factor, 50, probability=TRUE, xlab="N_UMI/<N_UMI>", main = main.use)
    } else {
        hist(size_factor, 50, probability=TRUE, xlab="N_Reads/<N_Reads>", main = main.use)
    }
    curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.999), add=TRUE, col="red", main="Gamma dist fit for size factor")
    text(5,0.6, paste("shape = ", round(fit$estimate[1],2)))
    text(5,0.5, paste("rate = ", round(fit$estimate[2],2)))
    }

    # Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
    # then cX ~ Gamma(a, b/c)
    a_i = rep(fit$estimate[1], length(mean_emp)); names(a_i) = names(mean_emp)
    b_i = fit$estimate[2] / mean_emp; names(b_i) = names(mean_emp)
    mean_NB = a_i / b_i; var_NB = a_i*(1+b_i) / (b_i^2)
    cv_NB = sqrt(var_NB)/mean_NB

    diffCV = log(cv_emp) - log(cv_NB)
    if(ret.diffCV == TRUE){
        return(diffCV)
    }
    print(length(diffCV))

    if(!is.null(diffCV.cutoff)){
	pass.cutoff=names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > 0.005 & mean_emp < 100))]
    } else if (!is.null(diffCV.num_genes)){
        pass.cutoff=names(sort(diffCV,decreasing=T)[1:diffCV.num_genes])
    }
    print(pass.cutoff)
    if(do.plot == T){
        plot(mean_emp,cv_emp,pch=16,cex=0.5,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy", main = main.use)
        if (do.spike) points(mean_emp[spike.genes],cv_emp[spike.genes],pch=16,cex=0.5,col="red")
        curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
        or = order(mean_NB)
        lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
        if(do.text) text(mean_emp[pass.cutoff],cv_emp[pass.cutoff],pass.cutoff,cex=cex.text.use)
        points(mean_emp[pass.cutoff], cv_emp[pass.cutoff], pch=16,cex=0.5,col="red")
        hist(diffCV, 50, probability=TRUE, xlab="diffCV")
        par(mfrow=c(1,1))
    }

    return(as.character(pass.cutoff))
}
