loess_regression = function(...){
    # fit loess curve and get fitted values
    fit = loess(...)
    p.x = fit$x[order(fit$x)]
    p.y = fit$fitted[order(fit$x)]
    return(list(fit=fit, x=p.x, y=p.y))
}


mixed_glm = function(x, y, k=2, ...){
    require(flexmix)
    # fit mixed glm and get assignments
    fit = flexmix(y ~ x, k=k, ...)
    groups = apply(fit@posterior[[1]], 1, which.max)
    return(list(fit=fit, groups=groups))
}


select_points = function(x, y, n, dir='both', loess=FALSE, nbins=25, bin_type='equal_width', exclude=c(0,0)){
    require(Hmisc)
    source('~/code/single_cell/downsample.r')

    # fix inputs
    if(!is.list(exclude)){exclude=list(exclude)}
    i = ((is.na(x) | is.na(y)) | (is.infinite(x) | is.infinite(y)))
    i = which(!i)
    xi = x[i]
    yi = y[i]

    # de-trend
    if(loess == TRUE){
        l = loess_regression(yi ~ xi, family='symmetric')
	yi = l$fit$residuals
    }

    # exclude data
    j = apply(sapply(exclude, function(e) (e[[1]] < xi) & (xi < e[[2]])), 1, any)
    j = which(!j)
    xi = xi[j]
    yi = yi[j]
    i = i[j]

    # bin x-axis
    if(bin_type == 'equal_width'){
        groups = cut2(xi, cuts=seq(from=min(xi, na.rm=T), to=max(xi, na.rm=T), length.out=nbins), m=2*n/nbins)
    } else {
        groups = cut2(xi, g=nbins)
    }

    # points
    j = c()

    # up points
    if(dir %in% c('up', 'both')){
        j = c(j, as.numeric(simple_downsample(cells=1:length(xi), groups=groups, ngene=yi, total_cells=n)))
    }

    # down points
    if(dir %in% c('down', 'both')){
        j = c(j, as.numeric(simple_downsample(cells=1:length(xi), groups=groups, ngene=-1*yi, total_cells=n)))
    }

    return(i[j])
}

