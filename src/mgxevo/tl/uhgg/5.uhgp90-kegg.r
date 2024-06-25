library(stringr)
library(data.table)

library(argparse)
parser = ArgumentParser()
parser$add_argument('--uhgp90_dir')
parser$add_argument('--output')
args = parser$parse_args()

out.dir  = args$output
uhgp.dir = args$uhgp90_dir

ifn = sprintf('%s/uhgp-90_eggNOG.tsv', uhgp.dir)

cmd = sprintf('grep "ko:" %s | cut -f 1,3,9', ifn)

df1 = fread(cmd = cmd,sep='\t',header=F)
colnames(df1) = c('id','id2','ko')

df1[,id2:= substr(id2, 1,22)]
df1 = melt(df1, id='ko')[,.(id=value, ko)]
df1[,ko := lapply(.SD, strsplit, split=','), .SDcols='ko']
df1 = df1[, lapply(.SD, unlist), by=1:nrow(df1)][,nrow:=NULL][]
df1[,ko:=substr(ko,4,11)]
df1 = unique(df1)

ofn = sprintf('%s/uhgp-90.ko', out.dir)
fwrite(df1,ofn,sep='\t',quote=F,row.names = F)