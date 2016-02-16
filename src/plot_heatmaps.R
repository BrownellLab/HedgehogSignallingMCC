## Heatmap of Hedgehog pathway
expr.data <- data.frame(exprs(gse.basalcell.hedgehog), ID=rownames(exprs(gse.basalcell.hedgehog)), gene.symbol=fData(gse.basalcell.hedgehog)[['Gene Symbol']])
expr.data <- transform(expr.data, plot.name=paste(gene.symbol, ID))
expr.data <- transform(expr.data, plot.name=factor(plot.name, levels=plot.name[order(as.character(gene.symbol))], ordered=TRUE))

my.mat.melted <- melt(expr.data, id.vars=c("ID", "gene.symbol", "plot.name"))
colnames(my.mat.melted) <- c("ID", "gene.symbol", "plot.name", "sample", "value")

expr.max <- round(max(exprs(gse.basalcell)) + 1)
expr.min <- round(min(exprs(gse.basalcell)) - 1)
expr.mid <- round(median(exprs(gse.basalcell)))

p <- ggplot(my.mat.melted, aes(y=plot.name, x=sample))
p <- p + geom_tile(aes(fill=value), size=0.5, color="white")
p <- p + scale_fill_gradient2(low="blue", mid="white", high="red", breaks=expr.min:expr.max, midpoint=expr.mid)
p <- p + labs(y="Probe", x="Sample")
p <- p + theme_bw()
p <- p + opts(axis.text.x=theme_text(size=11, angle=270), axis.text.y=theme_text(size=12)) ## + coord_flip()
## p <- p + facet_grid(cancertype + sample.type, scales="free", space="free")
## p <- p + facet_grid(. ~ geneclass, scales="free", space="free")

print(p)
ggsave("graphs/heatmap_basal.eps", width=10, height=10)
ggsave("graphs/heatmap_basal.png", width=10, height=10)
dev.off()

expr.data <- data.frame(exprs(gse.medullo.hedgehog), ID=rownames(exprs(gse.medullo.hedgehog)),
                        gene.symbol=fData(gse.medullo.hedgehog)[['Gene Symbol']])

expr.data <- transform(expr.data, plot.name=paste(gene.symbol, ID))
expr.data <- transform(expr.data, plot.name=factor(plot.name, levels=plot.name[order(as.character(gene.symbol))], ordered=TRUE))

my.mat.melted <- melt(expr.data, id.vars=c("ID", "gene.symbol", "plot.name"))
colnames(my.mat.melted) <- c("ID", "gene.symbol", "plot.name", "sample", "value")

# Add in subgroup
my.mat.melted <- merge(my.mat.melted,
                       with(pData(gse.medullo.hedgehog), data.frame(sample=sampleNames(gse.medullo.hedgehog), subgroup=gsub('subgroup: ', '', characteristics_ch1))),
                       by.x="sample", by.y="sample")

expr.max <- round(max(exprs(gse.medullo)) + 1)
expr.min <- round(min(exprs(gse.medullo)) - 1)
expr.mid <- round(median(exprs(gse.medullo)))

p <- ggplot(my.mat.melted, aes(y=plot.name, x=sample))
p <- p + geom_tile(aes(fill=value), size=0.5, color="white")
p <- p + scale_fill_gradient2(low="blue", mid="white", high="red", breaks=expr.min:expr.max, midpoint=expr.mid)
p <- p + labs(y="Probe", x="Sample")
p <- p + theme_bw()
p <- p + theme(axis.text.x=element_text(size=11, angle=270), axis.text.y=element_text(size=12)) ## + coord_flip()
## p <- p + facet_grid(cancertype + sample.type, scales="free", space="free")
p <- p + facet_grid(. ~ subgroup, scales="free", space="free")

print(p)
ggsave("graphs/heatmap_medullo.eps", width=20, height=10)
ggsave("graphs/heatmap_medullo.png", width=20, height=10)
dev.off()



expr.data <- data.frame(exprs(eset.hedgehog), ID=rownames(exprs(eset.hedgehog)), gene.symbol=fData(eset.hedgehog)[['gene.symbol']])
expr.data <- transform(expr.data, plot.name=paste(gene.symbol, ID))
expr.data <- transform(expr.data, plot.name=factor(plot.name, levels=plot.name[order(gene.symbol)], ordered=TRUE))

my.mat.melted <- melt(expr.data, id.vars=c("ID", "gene.symbol", "plot.name"))
colnames(my.mat.melted) <- c("ID", "gene.symbol", "plot.name", "sample", "value")

expr.max <- round(max(exprs(eset)) + 1)
expr.min <- round(min(exprs(eset)) - 1)
expr.mid <- round(median(exprs(eset)))

p <- ggplot(my.mat.melted, aes(y=plot.name, x=sample))
p <- p + geom_tile(aes(fill=value), size=0.5, color="white")
p <- p + scale_fill_gradient2(low="blue", mid="white", high="red", breaks=expr.min:expr.max, midpoint=expr.mid)
p <- p + labs(y="Probe", x="Sample")
p <- p + theme_bw()
p <- p + opts(axis.text.x=theme_text(size=11, angle=270), axis.text.y=theme_text(size=12)) ## + coord_flip()
## p <- p + facet_grid(cancertype + sample.type, scales="free", space="free")
## p <- p + facet_grid(. ~ geneclass, scales="free", space="free")

print(p)
ggsave("graphs/heatmap_mcc.eps", width=10, height=10)
ggsave("graphs/heatmap_mcc.png", width=10, height=10)
dev.off()



## p <- ggplot(my.mat.melted, aes(x=plot.name, y=value))
## p <- p + geom_boxplot(aes(fill=samplelevel))
## p <- p + facet_wrap(~ plot.name)

## p <- p + labs(x="Probe", y="Sample")
## p <- p + theme_bw()
## p <- p + opts(axis.text.x=theme_text(size=10, angle=270), axis.text.y=theme_text(size=12)) ## + coord_flip()
## p <- p + facet_grid(class + cellline ~ geneclass, scales="free", space="free")
## p <- p + facet_grid(samplelevel ~ geneclass, scales="free", space="free")
## ## p <- p + facet_grid(. ~ geneclass, scales="free", space="free")
## p <- p + coord_equal()

## ## pdf("graphs/heatmap_markers.pdf", width=10, height=5)
## print(p)
## ## dev.off()

