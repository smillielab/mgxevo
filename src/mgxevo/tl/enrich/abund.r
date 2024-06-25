script_path <- this.path::this.path()
script_dir <- dirname(script_path)

treetool <- paste0(script_dir,'/../tree/treetools.r')
util <- paste0(script_dir,'/../../ut/util.r')

source(treetool)
source(util)

library(argparse)

# create parser
parser <- ArgumentParser()

# add arguments
parser$add_argument("--bowtie_aln", help="Path to bowtie alignment file")
parser$add_argument("--name", help="Name of gene")
parser$add_argument("--metadata", help="Path to metadata file")
parser$add_argument("--multigene_table", help="Path to multigene table")
parser$add_argument("--output_dir", help="Path to output directory")
parser$add_argument("--min_gene_count", help="Minimum gene count for a given sample", type="integer", default=1000)
parser$add_argument("--min_sample_count", help="Minimum sample count for a given gene", type="integer", default=100)
parser$add_argument("--genes", type="character", nargs='+',
                        help='a list of marker genes to be processed')


# parse arguments
args <- parser$parse_args()
gn <- args$name
type <- args$type
genes <- args$genes

# ------------------------------
# load counts and other metadata
# ------------------------------

# load data
x = ffread(args$bowtie_aln, row.names=1, header=T, sep=',')
x = Matrix(as.matrix(x), sparse=T)

meta = read.metadata.nice(args$metadata)

# gene panel for abundance estimates
mx = ffread(args$multigene_table)
rownames(mx) = mx$genome

# filtering
# exclude samples with <1000 counts across all genes, and genes that occur in <100 samples
i = rowSums(x) >= args$min_gene_count
j = colSums(x > 0) >= args$min_sample_count
x = x[i,j]

# # fast row normalisation
x_row_names = rownames(x)
x = Matrix::Diagonal(x = 1 / Matrix::rowSums(x)) %*% x
x = log2(1e4*x + 1)
rownames(x) = x_row_names

# # aggregate by patient 
g = meta[match(rownames(x),meta$id),]$subject
x = sparse_agg(x, g, 'mean')

# # only look at genomes that have all marker genes present after filtering
mx = mx[apply(mx[,genes], 1, function(a) all(a %in% colnames(x))),]

# --------------
# run lmer model
# --------------


# IBD vs HC
ri = do.call(rbind, sapply(1:nrow(mx), function(i){
                df = rowMeans(x[,as.character(mx[i, genes])])
                df = cbind(counts = df, meta[match(names(df), meta$subject),c('cohort','nice_dx','nice_dx2','nice_age','nice_sex','nice_bmi')])
                df[,c('nice_age','nice_bmi')] = sapply(df[,c('nice_age','nice_bmi')] , unlist)
                df$nice_dx = factor(df$nice_dx, levels=c('HC', 'IBD'))
                df$nice_dx2 = factor(df$nice_dx2, levels=c('HC', 'UC', 'CD'))

                p = hush(lmer(counts ~ nice_dx + nice_age + nice_sex + nice_bmi + (1|cohort), data=df))
                u = summary(p)$coeff
                u = u[setdiff(rownames(u),'(Intercept)'),c('Estimate', 'Pr(>|t|)')]
                colnames(u) = c('coeff','pval')
                u = cbind(name = mx[i,]$genome, data.frame(u), test = rownames(u))
                return(u)},simplify=F))

# (UC|CD) vs HC
ri = rbind(ri, do.call(rbind, sapply(1:nrow(mx), function(i){
                df = rowMeans(x[,as.character(mx[i, genes])])
                df = cbind(counts = df, meta[match(names(df), meta$subject),c('cohort','nice_dx','nice_dx2','nice_age','nice_sex','nice_bmi')])
                df[,c('nice_age','nice_bmi')] = sapply(df[,c('nice_age','nice_bmi')] , unlist)
                df$nice_dx = factor(df$nice_dx, levels=c('HC', 'IBD'))
                df$nice_dx2 = factor(df$nice_dx2, levels=c('HC', 'UC', 'CD'))
                
                p = hush(lmer(counts ~ nice_dx2 + nice_age + nice_sex + nice_bmi + (1|cohort), data=df))
                u = summary(p)$coeff
                u = u[setdiff(rownames(u),'(Intercept)'),c('Estimate', 'Pr(>|t|)')]
                colnames(u) = c('coeff','pval')
                u = cbind(name = mx[i,]$genome, data.frame(u), test = rownames(u))
                u = u[u$test %in% c('nice_dx2CD', 'nice_dx2UC'),]
                return(u)},simplify=F)))

# CD vs other
ri = rbind(ri, do.call(rbind, sapply(1:nrow(mx), function(i){
                df = rowMeans(x[,as.character(mx[i, genes])])
                df = cbind(counts = df, meta[match(names(df), meta$subject),c('cohort','nice_dx','nice_dx2','nice_age','nice_sex','nice_bmi')])
                df[,c('nice_age','nice_bmi')] = sapply(df[,c('nice_age','nice_bmi')] , unlist)
                df$nice_dx = factor(df$nice_dx, levels=c('HC', 'IBD'))
                df$nice_dx2 = factor(df$nice_dx2, levels=c('HC', 'UC', 'CD'))
                df$other = df$nice_dx2
                df$other = factor(gsub('HC|UC','other', df$other), levels=c('other', 'CD'))

                p = hush(lmer(counts ~ other + nice_age + nice_sex + nice_bmi + (1|cohort), data=df))
                u = summary(p)$coeff
                u = u[setdiff(rownames(u),'(Intercept)'),c('Estimate', 'Pr(>|t|)')]
                colnames(u) = c('coeff','pval')
                u = cbind(name = mx[i,]$genome, data.frame(u), test = rownames(u))
                u = u[grepl('other', u$test),]
                return(u)},simplify=F)))

# UC vs other
ri = rbind(ri, do.call(rbind, sapply(1:nrow(mx), function(i){
                df = rowMeans(x[,as.character(mx[i, genes])])
                df = cbind(counts = df, meta[match(names(df), meta$subject),c('cohort','nice_dx','nice_dx2','nice_age','nice_sex','nice_bmi')])
                df[,c('nice_age','nice_bmi')] = sapply(df[,c('nice_age','nice_bmi')] , unlist)
                df$nice_dx = factor(df$nice_dx, levels=c('HC', 'IBD'))
                df$nice_dx2 = factor(df$nice_dx2, levels=c('HC', 'UC', 'CD'))
                df$other = df$nice_dx2
                df$other = factor(gsub('HC|CD','other', df$other), levels=c('other', 'UC'))

                p = hush(lmer(counts ~ other + nice_age + nice_sex + nice_bmi + (1|cohort), data=df))
                u = summary(p)$coeff
                u = u[setdiff(rownames(u),'(Intercept)'),c('Estimate', 'Pr(>|t|)')]
                colnames(u) = c('coeff','pval')
                u = cbind(name = mx[i,]$genome, data.frame(u), test = rownames(u))
                u = u[grepl('other', u$test),]
                return(u)},simplify=F)))

# ------------
# write output
# ------------
              
if(is.null(gn) & length(genes) > 1) {
    gn <- "multigene"
} else if (is.null(gn)) {
    gn <- genes
}
print(paste('Saving output as ', paste0(args$output_dir, '/', gn, '.abund.enrich.csv')))
fwrite(ri, paste0(args$output_dir, '/', gn, '.abund.enrich.csv'))