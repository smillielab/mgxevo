suppressMessages({
    treetools();
    singlecell()})

mx = ffread('~/smillie-data/db/uhgg/genes/ids.multigene.nice.txt')
nr = readLines('/broad/smillie-data/proj/mgxevo/meta/multigene.nr.txt')
rownames(mx) = mx$genome
ph = do.call(rbind, sapply(mx$genome, function(gi) fread(sprintf('/broad/smillie-data/proj/mgxevo/phylo/multigene/%s.csmillie.full.final.txt', gi)) ,simplify=F))
ph = ph[test == 'nice_dxIBD' | is.na(test)]
ph[, padj:=p.adjust(pval, 'fdr'), .(name)]
ph = ph[,.SD[which.min(padj)],.(name)]
ph[,genome:=name][,name:=NULL]
ph = merge(ph, mx, by='genome', all=T)
uf = ffread('/broad/smillie-data/proj/mgxevo/figures/multi.unifrac.tsv',as.dt=T)
ph = merge(ph, uf,  by='genome', all=T)
setkey(ph, genome)
ph = ph[nr]
                           
vv = ph[padj<.01,.(x=abs(coef),y=unifrac)]
pp = simple_scatter(vv$x,vv$y,
               xlab='Effect size',
               ylab='UniFrac',
               size=2) + 
    geom_smooth(aes(x=x,y=y),method = lm, se=FALSE, color='black') + 
    ggtitle(sprintf('Spearman\'s rho = %.2f',cor(vv$x, vv$y,use = 'pairwise.complete.obs',method = 'spearman')))
                           
pdf('unifrac_vs_glmer.pdf',5,5)
pp
dev.off()