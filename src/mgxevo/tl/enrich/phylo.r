script_path <- this.path::this.path()
script_dir <- dirname(script_path)

treetool <- paste0(script_dir,'/../tree/treetools.r')
util <- paste0(script_dir,'/../../ut/util.r')

source(treetool)
source(util)

library(lme4)
library(lmerTest)
library(glmnet)
library(argparse)

# create parser
parser <- ArgumentParser()

# add arguments
parser$add_argument("--tree_file", help="Path to tree file")
parser$add_argument("--name", help="Name of gene")
parser$add_argument("--type", help="type argument", default='full')
parser$add_argument("--metadata", help="Path to metadata file")
parser$add_argument("--output_dir", help="Path to output directory")

# parse arguments
args <- parser$parse_args()
gn <- args$name
type <- args$type

# build dataframe for regression
# ------------------------------

# process metadata
m = read.metadata.nice(args$metadata)

# # read tree
x = read.tree.nice(args$tree_file)
d = setRownames(m, m$id)
d = d[strip.dx_labs(x$tip.label),]

# run regression model on each node
# ---------------------------------

res = sapply(Ntip(x)+(1:Nnode(x)), function(node){
        
        # add node to input data
        d$in_node = d$id %in% strip.dx_labs(extract.clade(phy = x, node =  node)$tip.label)
	
	name = gn
	node = node
	test = NA
	coef = NA
	pval = NA
	note = NA
	
	# require minimum n=10 per node
	if((min(table(d$in_node)) > 10) && (min(table(d$nice_dx)) > 10)){
	
		tryCatch({
		    if(type == 'dx'){
		        u = summary(glmer(in_node ~ nice_dx + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))
			test = rownames(u$coef)
			coef = u$coef[,'Estimate']
			pval = u$coef[,'Pr(>|z|)']
			if(!is.null(u$optinfo$conv$lme4$messages)){note = u$optinfo$conv$lme4$messages}
		    }
		    if(type == 'full'){
			u = summary(glmer(in_node ~ nice_dx + nice_age + nice_sex + nice_bmi + (1|cohort), data=d, family=binomial(), glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=1e7))))
			test = rownames(u$coef)
			coef = u$coef[,'Estimate']
			pval = u$coef[,'Pr(>|z|)']
                        if(!is.null(u$optinfo$conv$lme4$messages)){note = u$optinfo$conv$lme4$messages}
		    }
	        },
	        error = function(...){
		    note = 'error'
		},
		warning=function(...){
		    note = 'warning'
		}
		)
	   	
	}
	ri = cbind(name=name, node=node, test=test, coef=coef, pval=pval, note=note)
	ri
	
}, simplify=F)

res = as.data.frame(do.call(rbind, res))
res$node = as.integer(res$node)
res$coef = as.numeric(res$coef)
res$pval = as.numeric(res$pval)
res = as.data.table(res)

write.table(res, paste0(args$output_dir, '/', gn,'.phylo.enrich.', type, '.txt'), row.names=F, quote=F, sep='\t')
