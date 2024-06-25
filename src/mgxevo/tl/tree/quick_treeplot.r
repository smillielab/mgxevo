script_path <- this.path::this.path()
script_dir <- dirname(script_path)

treetools <- paste0(script_dir,'/treetools.r')
source(treetools)

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--gn")
parser$add_argument("--tree_dir")
parser$add_argument("--enrich_dir")
args <- parser$parse_args()

gn <- args$gn
tree_dir <- args$tree_dir
enrich_dir <- args$enrich_dir

quick_treeplot(gn, tree_dir, enrich_dir)