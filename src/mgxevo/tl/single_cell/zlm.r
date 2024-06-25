load_mast()


# Fast way to run MAST using a dataframe and a formula
# Inputs:
# - data = data.frame containing all data (gene expression + covariates)
# - h1 = named list of hypotheses to test
#   example: list(MZB1 = NLRP7 ~ MZB1 + nGene)
# - h2 = named list of null hypotheses
#   example: list(gene = NLRP7 ~ nGene)
# Returns:
# - coefficients (coefD, coefC)
# - p-values (pvalD, pvalC, pvalH)


trycatch = function(a){tryCatch(expr={a}, error=function(e){NULL})}

run_zlm = function(data, h1, h0){

    # check hypotheses
    if(!is.list(h1)){h1 = list(h1)}
    if(!is.list(h0)){h0 = list(h0)}

    # track progress
    i = 0
    cat(paste(c('\nTesting ', length(h1), ' hypotheses\n', rep('-', 100), ' 100%', '\n'), collapse=''))

    # run zlm for every h1 hypothesis
    r1 = sapply(h1, function(f){
        i <<- i + 1
	if(as.integer(100*i/length(h1)) != as.integer(100*(i-1)/length(h1))){cat('.')}
        trycatch(zlm(f, sca=data))
    }, simplify=F)
    cat('\n\n')

    # run zlm for every h0 hypothesis
    r0 = sapply(h0, function(formula){
        trycatch(zlm(formula, sca=data))
    }, simplify=F)

    # get coefficients and p-values for every h1 hypothesis
    sapply(r1, function(a){

        coefD = trycatch(t(a$disc$coefficients))
	coefC = trycatch(t(a$cont$coefficients))
	pvals = do.call(cbind, sapply(r0, function(b) trycatch(zlm.pvals(a, b)), simplify=F))

	cnames = unique(c(colnames(coefD), colnames(coefC), colnames(pvals)))
	res = matrix(NA, nrow=5, ncol=length(cnames))
	rownames(res) = c('coefD', 'pvalD', 'coefC', 'pvalC', 'pvalH')
	colnames(res) = cnames

	res['coefD', colnames(coefD)] = coefD
	res['coefC', colnames(coefC)] = coefC
	res[rownames(pvals), colnames(pvals)] = pvals

	res
    }, simplify=F)
}

zlm.loglik = function(zlm_out){

    # Initialize log-likelihoods
    L = c(C=0, D=0)

    # Discrete
    dev = zlm_out$disc$deviance
    L['D'] = -dev/2

    # Continuous
    s2 = zlm_out$cont$dispersionMLE
    dev = zlm_out$cont$deviance
    N = zlm_out$cont$df.null+1
    L['C'] = -.5*N*(log(s2*2*pi)+1)

    return(L)
}

zlm.df = function(zlm_out){
    c(C=zlm_out$cont$df.residual, D=zlm_out$disc$df.residual)
}

zlm.pvals = function(r1, r0){

    # log-likelihoods
    l0 = zlm.loglik(r0)
    l1 = zlm.loglik(r1)

    # df.residuals
    d0 = zlm.df(r0)
    d1 = zlm.df(r1)

    # p-value calculation
    lambda = -2*(l0 - l1)
    df = d0 - d1
    pvals = MAST:::makeChiSqTable(lambda, df, 1)[,'Pr(>Chisq)']
    names(pvals) = list(cont='pvalC', disc='pvalD', hurdle='pvalH')[names(pvals)]

    return(pvals)
}
