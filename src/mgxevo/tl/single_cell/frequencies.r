
dirichlet_regression = function(counts, covariates, formula){

    # Dirichlet multinomial regression to detect changes in cell frequencies
    # formula is not quoted, example: counts ~ condition
    # counts is a [samples x cell types] matrix
    # covariates holds additional data to use in the regression
    #
    # Example:
    # counts = do.call(cbind, tapply(obj$meta.data$orig.ident, obj$ident, table))
    # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
    # res = dirichlet_regression(counts, covariates, counts ~ condition)

    require(DirichletReg)

    # Fix counts matrix if generated with table
    counts = as.data.frame.matrix(counts)

    # Calculate regression
    counts$counts = DR_data(counts)
    data = cbind(counts, covariates)
    fit = DirichReg(counts ~ condition, data)

    # Get p-values
    u = summary(fit)
    pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
    v = names(pvals)
    pvals = matrix(pvals, ncol=length(u$varnames))
    rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
    colnames(pvals) = u$varnames
    fit$pvals = pvals

    return(fit)
}


