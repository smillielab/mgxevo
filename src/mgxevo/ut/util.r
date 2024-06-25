ffread = function(a, sep='\t', row.names=NULL, header=TRUE, as.dt=FALSE){
    library(data.table)
    if(grepl('gz', a)){
        h = strsplit(readLines(a, n=1), sep)[[1]]
	x = fread(paste0('zcat ', a), sep=sep, skip=1, header=F)
    } else {
        h = strsplit(readLines(a, n=1), sep)[[1]]
	x = fread(a, skip=1, sep=sep, header=F)
    }
    if(length(h) == (ncol(x)-1)){
        print('Setting header=T and row.names = 1')
	header = TRUE
        row.names = 1
    }
    # print(paste('ffread: reading', a, 'with header =', header, 'and row.names=', row.names))
    if(nrow(x) == 1){
        if(row.names == 1){
	    x = read.table(a, sep=sep, header=F, row.names=row.names, skip=1)
	}
    } else {
        x = data.frame(x, row.names=row.names)
    }
    if(length(h) == (ncol(x) + 1)){
        h = h[2:length(h)]
    }
    if(header == FALSE){
        if(length(h) == ncol(x)){
            x = rbind(h, x)
	} else {
	    stop('error: file dimensions')
	}
    } else {
        if(length(h) == ncol(x)){
            colnames(x) = h
        } else if(length(h) == ncol(x)+1){
            colnames(x) = h[2:length(h)]
        } else {
            stop('error: file dimensions')
        }
    }
    if(as.dt == TRUE){x = as.data.table(x)}
    return(x)
}

mapv = function(a,b){
    x = unique(data.frame(a=a, b=b, stringsAsFactors=F))
    setNames(x$b, x$a)
}

read.metadata.nice = function(ifn){
    # process metadata
    m = ffread(ifn)
    fix.names = c('Study cohort' = 'cohort',
                     'Subject' = 'subject',
                     'Sample ID' = 'sample',
                     'Diagnosis' = 'dx',
                     'Age (years)' = 'age',
                     'Sex' = 'sex',
                     'BMI (kg/m^2)' = 'bmi',
                     'Is metagenome new?' = 'not_published',
                     'Fecal calprotectin' = 'fcal',
                     'File path' = 'file',
                     'Paired?' = 'paired',
                     'Diagnosis (UC, CD, IBDU, HC)' = 'nice_dx',
                     'Deidentified file path' = 'file.renamed',
                     'Time relative to first sample (day; if longitudinal)' = 'day')
    nice2fancy = setNames(names(fix.names),(fix.names))
    colnames(m) = (fix.names[colnames(m)])
    m$id = m$sample
    m$nice_dx = factor(gsub("IBDU", "CD",m$nice_dx), levels=c('HC', 'CD', 'UC'))
    m$nice_dx2 = m$nice_dx
    m$nice_dx = factor(gsub("CD|UC", "IBD",m$nice_dx), levels=c('HC', 'IBD'))
    
    m$sex[m$sex == ''] = NA
    m$age[grepl('[A-Z]', m$age)] = NA
    m$age[m$age == ''] = NA
    
    m$age = as.numeric(m$age)
    m$sex = factor(m$sex, levels=c('M', 'F'))

    # impute missing values
    if(sum(!is.na(m$age)) < 10){m$age = rnorm(mean=0, sd=1, nrow(m))}
    if(sum(!is.na(m$bmi)) < 10){m$bmi = rnorm(mean=0, sd=1, nrow(m))}
    if(sum(!is.na(m$sex)) < 10){m$sex = sample(c('M', 'F'), nrow(m), replace=T)}
    m$nice_age = scale(ifelse(is.na(m$age), mean(m$age, na.rm=T), m$age))
    m$nice_bmi = scale(ifelse(is.na(m$bmi), mean(m$bmi, na.rm=T), m$bmi))
    m$nice_sex = m$sex
    i = is.na(m$nice_sex)
    m$nice_sex[i] = sample(m$nice_sex[!i], sum(i), replace=T)
    m$nice_sex = factor(m$nice_sex, levels=c('M', 'F'))
    return(m)
}



calc_rocr = function(predictions, labels, measures='auc', retx=FALSE){
    require(ROCR)

    if(length(measures) == 1){
        q = performance(prediction(predictions, labels, label.ordering=levels(labels)), measures)
    } else {
        q = performance(prediction(predictions, labels, label.ordering=levels(labels)), measures[[1]], measures[[2]])
    }

    x = q@x.values
    if(length(x) > 0){x = ifelse(is.na(x[[1]]), 0, x[[1]])}
    y = ifelse(is.na(q@y.values[[1]]), 0, q@y.values[[1]])

    if(length(y) > 1){
        if(q@x.name == 'Cutoff'){
	    i = which.max(y[x != min(x)])
	    x = x[i]
	    y = y[i]
	} else {
            y = sum(diff(x) * (head(y,-1) + tail(y,-1)))/2
	}
    }

    if(retx) {c(x,y)} else {round(y, 3)}
}


de.rocr = function(data, labels, measures='auc'){

    # Calculate AUC of ROC curve and average difference
    scores = sapply(rownames(data), function(a){calc_rocr(as.numeric(data[a,]), labels, measures=measures)})

    # Return marker genes
    markers = data.frame(gene=rownames(data), stringsAsFactors=FALSE)
    cname = paste(measures, collapse='_')
    markers[,cname] = scores
    markers = markers[order(markers[,cname], decreasing=T),]
    return(markers)
}


calc_fdr = function(predictions, labels){

    # Get classification cutoff with F measure
    f = calc_rocr(as.numeric(predictions), labels, 'f', retx=T)

    # Get true/false negative/positives
    u = factor(as.numeric(predictions >= f[[1]]), levels=c(0,1), ordered=T)
    q = as.vector(table(u, labels))
    tn = q[[1]]; fp = q[[2]]; fn = q[[3]]; tp = q[[4]];

    # Calculate statistics
    fdr = fp/(tp+fp)
    tpr = tp/(tp+fn)
    tnr = tn/(tn+fp)
    acc = (tp + tn)/sum(q)
    return(c(f[[1]], acc, tpr, tnr, fdr, f[[2]]))
}


de.fdr = function(data, labels, sens_cut=.1){

    # Calculate classification statistics
    markers = t(sapply(rownames(data), function(a){calc_fdr(data[a,], labels)}))
    colnames(markers) = c('cutoff', 'accuracy', 'sensitivity', 'specificity', 'fdr', 'f1')

    # Calculate average difference
    avg_diff = sapply(rownames(data), function(a){as.numeric(diff(tapply(as.numeric(data[a,]), labels, mean)))})

    # Return marker genes
    markers = cbind(markers, data.frame(avg_diff=avg_diff, gene=rownames(markers), stringsAsFactors=F))
    markers = markers[markers$sens >= sens_cut,]
    return(markers)

}

