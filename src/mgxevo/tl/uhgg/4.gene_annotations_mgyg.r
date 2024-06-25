script_path <- this.path::this.path()
script_dir <- dirname(script_path)

library(stringr)
library(data.table)

library(argparse)
parser = ArgumentParser()
parser$add_argument("--MGYG_id")
parser$add_argument('--uhgg_catalogue')
args = parser$parse_args()
mg = args$MGYG_id
uhgg_catalogue = args$uhgg_catalogue


ifn = paste0(uhgg_catalogue,str_sub(mg,1,13),'/',mg,'/genome/',mg,'.gtf')
ofn = paste0(uhgg_catalogue,str_sub(mg,1,13),'/',mg,'/genome/',mg,'.anno.nice.txt')

fl = readLines(ifn)
map = data.frame(do.call(rbind,sapply(fl, function(line){
    line = str_split(line,'\t')[[1]]
    if(line[3] %in% c('exon','gene')){
        return(c(refseq=str_match(line[9], 'gene_id "(.*?)"')[,2], gene=str_match(line[9], 'locus_tag "(.*?)"')[,2]))
    } else {
        return(c(refseq=NA, gene=NA))
    }}, simplify=F)))
rownames(map) = NULL
map = map[!duplicated(map[,1]),]
map = map[complete.cases(map),]


# load EggNOG

df1 = fread(cmd=paste0('python ', script_dir, '/tools/clean_eggNOG.py ' , uhgg_catalogue,str_sub(mg,1,13),'/',mg,'/genome/',mg,'_eggNOG.tsv'))

# tidy InterPro annotations
tidy.interpro = function(df2){
    df2 = unique(df2)
    gx = unique(df2$gene)
    do.call(rbind,sapply(gx, function(A) apply(df2[gene==A],2, 
                                           function(B){
                                               a = B[B!='']
                                               paste(unique(a),collapse=',')}),simplify=F))}

df2 = fread(cmd=paste0('python ', script_dir, '/tools/clean_interPro.py ', uhgg_catalogue, str_sub(mg,1,13),'/',mg,'/genome/',mg,'_InterProScan.tsv'))
df2 = tidy.interpro(df2)
                         
dfx = merge(df1, df2, by = 'gene', all.x=TRUE, all.y=TRUE)
dfx = do.call(rbind, sapply(unique(dfx$gene), function(gx){
    A = dfx[gene==gx]
    go.1 = stringr::str_split(A$go,',')[[1]]
    go.1 = unique(go.1[go.1!=''])
    
    go.2 = stringr::str_split(A$GO,',')[[1]]
    go.2 = unique(go.2[go.2!=''])
    
    go = c(go.1, go.2)
    go = unique(go[go!=''])
    go = paste(go, collapse=',')
    
    pathway.1 = stringr::str_split(A$pathway,',')[[1]]
    pathway.1 = unique(pathway.1[pathway.1!=''])
    
    pathway.2 = stringr::str_split(A$kegg_pathway,',')[[1]]
    pathway.2 = unique(pathway.2[pathway.2!=''])
    pathway.2 = paste0('KEGG:',pathway.2)
    
    pathway = c(pathway.1, pathway.2)
    pathway = unique(pathway[pathway!=''])
    pathway = setdiff(pathway,c('NA','KEGG:'))
    pathway = paste(pathway, collapse=', ')

    A$go = go
    A$pathway = pathway
    A$GO = NULL
    A$eval = 
    return(A)
},simplify = F))
res = merge(map,dfx, by='gene', all.x=TRUE)

write.table(res, file=ofn, row.names=F, quote=F, sep='\t')