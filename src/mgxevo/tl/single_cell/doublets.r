score_doublets = function(obj, sig, scores=NULL, doublets=NULL, pd=0.25, min_exp=1, min_odds=1, do.fast=FALSE){

    # fix input arguments
    if(is.null(doublets)){doublets = rep(FALSE, ncol(obj$data))}
    doublets[is.na(doublets)] = FALSE
    print(paste('Filtering', sum(doublets), 'doublets'))
    print(paste('Ignoring signatures with log2(TP10K+1) <', min_exp))
    print(paste('Odds ratio cutoff =', min_odds))

    # check column names
    sig = sig[names(sig) %in% levels(obj$ident)]
    if(length(sig) == 0){print('error: names(sig) != levels(obj$ident)')}

    # score cells
    if(is.null(scores)){
        print('Calculating signature scores')
	scores = score_cells(obj, names=sig, make.names=T)
    }
    i = !doublets
    scores.u = data.frame(aggregate(scores[i,], list(make.names(obj$ident[i])), mean), row.names=1)
    scores.s = data.frame(aggregate(scores[i,], list(make.names(obj$ident[i])), sd), row.names=1)

    # parameters for null hypothesis (h0)
    h0.u = scores.u[make.names(as.character(obj$ident)),]
    h0.s = scores.s[make.names(as.character(obj$ident)),]

    # parameters for alternative hypothesis (h1)
    h1.u = t(matrix(diag(as.matrix(scores.u[colnames(scores.u),])), nrow=ncol(scores.u), ncol=nrow(scores)))
    h1.s = t(matrix(diag(as.matrix(scores.s[colnames(scores.s),])), nrow=ncol(scores.s), ncol=nrow(scores)))
    h1.u = .5*h0.u + .5*h1.u
    h1.s = .5*sqrt(h0.s**2 + h1.s**2)

    if(do.fast == FALSE){

        # probabilities for null hypothesis (h0)
	print('Calculating H0 probabilities')
    	p0 = sapply(1:ncol(scores), function(j) {
      	    print(paste('Signature', colnames(scores)[j]))
    	    sapply(1:nrow(scores), function(i) {
                dnorm(scores[i,j], mean=h0.u[i,j], sd=h0.s[i,j])
    	    })
	})

	# probabilities for alternative hypothesis (h1)
	print('Calculating H1 probabilities')
	p1 = sapply(1:ncol(scores), function(j) {
    	    print(paste('Signature', colnames(scores)[j]))
    	    sapply(1:nrow(scores), function(i) {
                dnorm(scores[i,j], mean=h1.u[i,j], sd=h1.s[i,j])
    	    })
	})

    } else {

        library(dnormpar2)

	# probabilities for null hypothesis (h0)
	print('Fast calculating H0 probabilities')
	p0 = sapply(1:ncol(scores), function(j) {
	    print(paste('Signature', colnames(scores)[j]))
	    dnormpar2(scores[,j], h0.u[,j], h0.s[,j])
	})

	# probabilities for alternative hypothesis (h1)
	print('Fast calculating H1 probabilities')
	p1 = sapply(1:ncol(scores), function(j) {
	    print(paste('Signature', colnames(scores)[j]))
	    dnormpar2(scores[,j], h1.u[,j], h1.s[,j])
	})
    }

    # fix names
    rownames(h0.u) = rownames(h0.s) = rownames(h1.u) = rownames(h1.s) = rownames(p0) = rownames(p1) = rownames(scores)
    colnames(h0.u) = colnames(h0.s) = colnames(h1.u) = colnames(h1.s) = colnames(p0) = colnames(p1) = colnames(scores)

    # calculate odds ratio
    print('Detecting doublets')
    odds = (p1*pd)/(p0*(1-pd))
    odds[scores < min_exp] = 0

    # infer doublets
    doublets = apply(odds, 1, max) > min_odds

    # return data
    list(scores=scores, scores.u=scores.u, scores.s=scores.s, h0.u=h0.u, h0.s=h0.s, h1.u=h1.u, h1.s=h1.s, p0=p0, p1=p1, odds=odds, doublets=doublets)

}

find_doublets = function(obj, sig, scores=NULL, doublets=NULL, pd=0.25, min_exp=1, min_odds=c(4,2,1,1), do.fast=FALSE){

    print('Calculating signature scores')
    sig = sig[names(sig) %in% levels(obj$ident)]
    if(is.null(scores)){scores = score_cells(obj, names=sig)}
    if(is.null(doublets)){doublets = rep(FALSE, ncol(obj$data))}
    scores = scores[,make.names(names(sig))]

    res = list()
    for(i in 1:length(min_odds)){
        res[[i]] = score_doublets(obj=obj, sig=sig, scores=scores, doublets=doublets, pd=pd, min_exp=min_exp, min_odds=min_odds[[i]], do.fast=TRUE)
        doublets = res[[i]]$doublets
    }
    return(res)
}

#scores = score_cells(obj, sig=sig.use)
#res = find_doublets(obj, sig.use, scores=scores, pd=.25, min_exp=1.5, min_odds=c(8,7,6,5,4,3,2,1,4))
