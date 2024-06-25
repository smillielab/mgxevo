

# load_maps = function(h='hg19', m='mm10'){

#     # Get gene lists
#     hgenes = readLines(paste0('~/code/db/', h, '_genes.txt'))
#     mgenes = readLines(paste0('~/code/db/', m, '_genes.txt'))

#     # Get gene synonyms
#     hsyn = paste0('~/code/db/', h, '.gene_map.txt')
#     hsyn = read.table(hsyn, row.names=1, sep='\t', stringsAsFactors=F, quote='', comment.char='')
#     msyn = paste0('~/code/db/', m, '.gene_map.txt')
#     msyn = read.table(msyn, row.names=1, sep='\t', stringsAsFactors=F, quote='', comment.char='')

#     # Get other symbols
#     hsym = paste0('~/code/db/human.db2sym.txt')
#     hsym = read.table(hsym, row.names=1, sep='\t', stringsAsFactors=F, quote='', comment.char='')
#     msym = paste0('~/code/db/mouse.db2sym.txt')
#     msym = read.table(msym, row.names=1, sep='\t', stringsAsFactors=F, quote='', comment.char='')

#     # Combine synonyms and symbols
#     hsyn = unique(rbind(hsyn, hsym))
#     msyn = unique(rbind(msyn, msym))

#     # Get orthologs
#     h2m = paste0('~/code/db/orthologs.', h, '_to_', m, '.txt')
#     h2m = read.table(h2m, sep='\t', stringsAsFactors=F, quote='', comment.char='', row.names=1)
#     m2h = paste0('~/code/db/orthologs.', m, '_to_', h, '.txt')
#     m2h = read.table(m2h, sep='\t', stringsAsFactors=F, quote='', comment.char='', row.names=1)

#     # Return maps
#     return(list(h=h, m=m, hgenes=hgenes, mgenes=mgenes, hsyn=hsyn, msyn=msyn, hsym=hsym, msym=msym, h2m=h2m, m2h=m2h))
# }


# maps = load_maps()
# list2env(maps, .GlobalEnv)


predict_organism = function(genes){
    if(sum(genes %in% hgenes) > sum(genes %in% mgenes)){
        return('human')
    }
    if(sum(genes %in% hgenes) < sum(genes %in% mgenes)){
        return('mouse')
    }
    return('unknown')
}


fix_names = function(names){
    names = toupper(names)
    names = gsub('[^a-zA-Z0-9]', '', names)
    return(names)
}


get_synonyms = function(genes, target='human', do.unlist=TRUE){
    genes = fix_names(genes)
    if(target == 'human'){genes = hsyn[genes,1]}
    if(target == 'mouse'){genes = msyn[genes,1]}
    if(do.unlist == TRUE){
        genes = unlist(strsplit(genes, ','))
    } else {
        genes = strsplit(genes, ',')
    }
    genes
}

get_orthologs = function(genes, source='mouse', target='human'){
    genes = setNames(get_synonyms(genes, target=source, do.unlist=F), genes)
    genes = lapply(genes, fix_names)
    if(target == 'mouse'){ortho = lapply(genes, function(a) h2m[a,1])}
    if(target == 'human'){ortho = lapply(genes, function(a) m2h[a,1])}
    ortho = lapply(ortho, get_synonyms, target=target, do.unlist=T)
    ortho = lapply(ortho, function(a) sort(unique(a)))
    ortho
}


map_gene = function(genes, target='human', source='auto', do.unlist=TRUE){

    # predict source organism from gene list
    # --------------------------------------

    if(source == 'auto'){
        source = predict_organism(genes)
    }

    if(source == 'unknown'){
        cat(paste('\nmap_gene: unknown source organism, assuming source =', target, '\n'))
	source = target
    }

    # map genes using synonyms or orthologs
    # -------------------------------------

    if(source == target){
        get_synonyms(genes, target=target, do.unlist=do.unlist)
    } else {
        get_orthologs(genes, source=source, target=target)
    }

}


get_functions = function(genes, type='desc', org='auto'){
    library(rentrez)
    library(biomaRt)
    res = setNames(rep('', length(genes)), genes)
    if(org == 'auto'){
        org = predict_organism(genes)
    }
    if(type == 'desc'){
        if(org == 'human'){
            hmart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	    anno = getBM(attributes=c('hgnc_symbol', 'description'), mart=hmart)
    	} else if(org == 'mouse'){
            mmart = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
	    anno = getBM(attributes=c('mgi_symbol', 'description'), mart=mmart)
    	} else {
            stop('unknown organism')
    	}
    	anno = anno[anno[,1] != '',]
    	anno = anno[!duplicated(anno[,1]),]
    	anno = data.frame(anno, row.names=1)
    	anno[,1] = gsub('\\[.*\\]$', '', anno[,1])
	i = intersect(rownames(anno), genes)
	res[i] = anno[i,1]
    }
    if(type == 'summary'){
        library(mygene)
        if(org == 'human'){
	    library(org.Hs.eg.db)
	    entrez = mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
	    res[] = unlist(getGenes(entrez, fields='summary')$summary)
	} else if(org == 'mouse'){
	    library(org.Mm.eg.db)
	    entrez = mapIds(org.Mm.eg.db, genes, 'ENTREZID', 'SYMBOL')
	    res[] = unlist(getGenes(entrez, fields='summary')$summary)
	} else {
	    stop('unknown organism')
	}
    }
    return(res)
}
