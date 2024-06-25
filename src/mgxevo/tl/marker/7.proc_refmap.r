script_path <- this.path::this.path()
script_dir <- dirname(script_path)

library(data.table)
library(argparse)

util <- file.path(script_dir, "/../../ut/util.r")
source(util)

# marker panel 
library(argparse)

parser <- ArgumentParser(description = "Process input files.")
parser$add_argument("--multigene_table", help = "Path to multigene table.")
parser$add_argument("--blast_results", help = "Path to blast results file.")
parser$add_argument("--marker_mapping", help = "Path to marker mapping file.")
parser$add_argument("--output_dir", help = "Path to output directory.")
args <- parser$parse_args()

marker = fread(args$multigene_table)

# blast results
b6 = fread(args$blast_results)
colnames(b6) = c('query', 'target', 'id', 'alnlen', 'mism', 'opens', 'qlo', 'qhi', 'tlo', 'thi', 'evalue', 'bits')

m = fread(args$marker_mapping,header=T)
m = data.frame(m)
rownames(m) = m$genome
m$genome = NULL

# remove empty
i = apply(m != '', 1, all)
m = m[i,]

# split
i = rownames(m)
m = sapply(m, strsplit, ',')
rownames(m) = i


# iterate through combinations that are from the same genome
# NOTE: gsub('\\_NODE(.*)|^blast/','',a), needs to be replaced with the appropriate REGEX to get the orignal genome
v = sapply(rownames(m), function(ri){
        mi = m[ri,]
        index = Reduce(intersect, sapply(mi, function(a) unlist(gsub('_[^_]*$','',a)),simplify=F))
        u = t(sapply(index, function(id) sapply(names(mi), function(mx){
                ret = grep(id,mi[[mx]],value=T)
                if(length(ret)>1){
                    qry = unlist(marker[marker$genome==ri,mx,with=F])
                    hit = ret
                    df = b6[(b6$query == qry) & (b6$target %in% hit),]
                    ret = df[which.max(df$id),'target']
                }
                return(ret)})))
        u = cbind(genome = ri, u)
        rownames(u) = NULL
        return(u)
}, simplify = F)

# remove genomes where we don't have hits from the same org
i = sapply(v, dim)[2,] != 1
vq = v[i]
vq = do.call(rbind,vq)
fwrite(vq, paste0(args$output_dir, '/tree-seeds.mgyg.filtered.ref_map.proc.txt'), quote = FALSE, sep='\t', row.names=FALSE, col.names=FALSE)