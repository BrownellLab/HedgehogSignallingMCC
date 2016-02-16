## Merge datasets from hgu133plus2 and hgu133av2 by combining probes in common
common.probes <- intersect(rownames(fData(eset)), rownames(fData(gse.basalcell)))

## MCC data
eset.common <- eset[common.probes, ]

gse.basalcell.common <- gse.basalcell[common.probes, ]
gse.medullo.common <- gse.medullo[common.probes, ]
gse.sclc.common <- gse.sclc[common.probes, ]

samplenames <- c(sampleNames(eset), sampleNames(gse.medullo.common), sampleNames(gse.basalcell.common), sampleNames(gse.sclc.common))
pdata <- data.frame(sample=samplenames,
                    cancertype=c(as.character(pData(eset)$cancertype), rep("medullo", length(sampleNames(gse.medullo.common))), rep("basal", length(sampleNames(gse.basalcell.common))), rep("SCLC", length(sampleNames(gse.sclc.common)))),
                    batch=c(pData(eset)$batch, rep(3, length(sampleNames(gse.medullo.common))), rep(4, length(sampleNames(gse.basalcell.common))), rep(5, length(sampleNames(gse.sclc.common)))))

rownames(pdata) <- samplenames

combined.exprs <- cbind(exprs(eset), exprs(gse.medullo.common), exprs(gse.basalcell.common), exprs(gse.sclc.common))

neweset <- ExpressionSet(assayData=combined.exprs, phenoData=AnnotatedDataFrame(pdata))

## SLOW!!!!
neweset.exprs.batcheffect <- ComBat(dat=exprs(neweset), batch=pData(neweset)$batch, mod=NULL, par.prior=FALSE)

neweset.batcheffect <- ExpressionSet(assayData=neweset.exprs.batcheffect, phenoData=AnnotatedDataFrame(pData(neweset)))


## Z-score normalization
eset.common.exprs.zscore <- apply(exprs(eset.common), 1, function(x) (x - mean(x)) / sd(x))
gse.basalcell.common.exprs.zscore <- apply(exprs(gse.basalcell.common), 1, function(x) (x - mean(x)) / sd(x))
gse.medullo.common.exprs.zscore <- apply(exprs(gse.medullo.common), 1, function(x) (x - mean(x)) / sd(x))
gse.sclc.common.exprs.zscore <- apply(exprs(gse.sclc.common), 1, function(x) (x - mean(x)) / sd(x))

all.exprs.zscore <- t(rbind(eset.common.exprs.zscore, gse.medullo.common.exprs.zscore, gse.basalcell.common.exprs.zscore, gse.sclc.common.exprs.zscore))
neweset.common.zscores <- ExpressionSet(assayData=all.exprs.zscore, phenoData=AnnotatedDataFrame(pdata))

neweset.common.zscores <- ExpressionSet(assayData=t(eset.common.exprs.zscore), phenoData=AnnotatedDataFrame(pData(eset.common)))
