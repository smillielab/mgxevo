# input:
#   mi = MGYG id
# requires:
#   uhgg.sqlite: generated using the following SQLite commands:
#        sqlite3 uhgg.sqlite
#        sqlite> create table map(seed,gene);
#        sqlite> .mode tabs
#        sqlite> .import /broad/smillie-data/db/uhgg/uhgp-90/uhgp-90.tsv map
#        sqlite> .schema map
# output:
#   *.locus.mapped.csv: [genes x genomes] = UHGP genes
#   *.uhgp.csv = UHGP genes (for MGYG IDs that do not have a pan-genome)

library(dplyr)
library(dbplyr)
library(data.table)
library(stringr)

args = commandArgs(trailingOnly=T)
mi = args[[1]]
uhgg_path = args[[2]]

# get filepath and load matrix IFF it exists
ifn = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/genes_presence-absence_locus.csv', uhgg_path, substr(mi,1,13), mi)

if(file.exists(ifn)){
    y = fread(ifn)
    
    # what genes are in the matrix
    j = grep('GUT_GENOME\\d+', colnames(y), value=T)
    yq = data.frame(y[,..j])
    has.split = sapply(1:ncol(yq), function(j) grepl('\t', yq[,j]))
    genes = yq[!has.split]
    tmp = yq[has.split]
    tmp = str_split(tmp,"\t")
    tmp = unlist(tmp, use.names=F)
    genes = c(genes, tmp)
    genes = unique(genes)
    genes = setdiff(genes, c('',NA, 'NA'))
    
    # SQL query
    db = DBI::dbConnect(RSQLite::SQLite(), sprintf("%s/uhgg.sqlite", uhgg_path))
    qry = paste(sprintf('"%s"',unique(genes)), collapse=',')
    qry = sprintf("SELECT * FROM map WHERE gene IN (%s)", qry)
    res = collect(tbl(db, sql(qry)))
    
    # map genes onto seed values
    mp = setNames(res$seed, res$gene)
    
    zq = data.frame(y[,-c(1:14)])
    has.split = sapply(1:ncol(zq), function(j) grepl('\t', zq[,j]))
    zq[!has.split] = mp[zq[!has.split]]
    if(any(has.split)){
        tmp = zq[has.split]
        tmp = str_split(tmp,"\t")
        max.length = max(sapply(tmp, length))
        tmp = lapply(tmp, function(v) { c(v, rep(NA, max.length-length(v)))})
        tmp = matrix(unlist(tmp, use.names = FALSE), ncol = max.length, byrow = TRUE)
        tmp = do.call(rbind, sapply(1:ncol(tmp), function(j) mp[tmp[,j]], simplify=FALSE))
        tmp = sapply(1:ncol(tmp), function(j){
                    z = unique(tmp[,j])
                    z = z[!is.na(z)]
                    z = paste(z, collapse="\t")
                    return(z)
                })
        zq[has.split] = tmp
    }
    yq = cbind(y[,1:14], zq)
                                
    ofn = sprintf('%s/uhgg_catalogue/%s/%s/pan-genome/genes_presence-absence_locus.mapped.csv', uhgg_path, substr(mi,1,13), mi)
    write.table(yq, ofn, sep=",", row.names = F)
} else {
    # get genes
    ifn = sprintf('%s/uhgg_catalogue/%s/%s/genome/%s.gff', uhgg_path, substr(mi,1,13), mi, mi)
    cmd = sprintf("cut -d'=' -f2  %s | cut -d';' -f1 | grep -P '_[0-9]+$'", ifn)
    genes = unique(system(cmd, intern=T))

    # SQL query
    db = DBI::dbConnect(RSQLite::SQLite(), sprintf("%s/uhgg.sqlite", uhgg_path))
    qry = paste(sprintf('"%s"',unique(genes)), collapse=',')
    qry = sprintf("SELECT * FROM map WHERE gene IN (%s)", qry)
    res = collect(tbl(db, sql(qry)))

    # write files
    ofn = sprintf('%s/uhgg_catalogue/%s/%s/genome/%s.uhgp.txt', uhgg_path, substr(mi,1,13), mi, mi)
    writeLines(unique(res$seed), ofn)
}