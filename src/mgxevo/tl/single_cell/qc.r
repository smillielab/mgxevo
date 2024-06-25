

mem_apply = function(data, axis=2, FUN, k=25){
    library(Hmisc)
    x = c()
    i = 0
    if(axis == 1){
        g = cut2(1:nrow(data), g=k)
	for(gi in levels(g)){cat(paste0(100*i/k, '%'), ''); i = i + 1;
	    x = c(x, apply(data[g == gi,], 1, FUN))
	}
    } else {
        g = cut2(1:ncol(data), g=k)
	for(gi in levels(g)){cat(paste0(100*i/k, '%'), ''); i = i + 1;
	    x = c(x, apply(data[,g == gi], 2, FUN))
	}
    }
    cat('100%')
    x
}

calc_entropy = function(obj, groups=25){
    mem_apply(data=obj$data, axis=2, FUN=function(a){a=a/sum(a); sum(-a*log(a), na.rm=T)}, k=groups)
}
