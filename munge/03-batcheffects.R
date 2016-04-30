combined.exprs <- cbind(exprs(eset), exprs(gse.medullo.common), exprs(gse.basalcell.common), exprs(gse.sclc.common))
neweset <- ExpressionSet(assayData=combined.exprs, phenoData=AnnotatedDataFrame(pdata))

## SLOW!!!!
neweset.exprs.batcheffect <- ComBat(dat=exprs(neweset), batch=pData(neweset)$batch, mod=NULL, par.prior=FALSE)

neweset.batcheffect <- ExpressionSet(assayData=neweset.exprs.batcheffect, phenoData=AnnotatedDataFrame(pData(neweset)))

## Do quality control on batch corrected data
qc.all.batcheffect <- arrayQualityMetrics(expressionset=neweset.batcheffect, outdir="reports/QC/RMA.batcheffect.all/", force=T, intgroup="cancertype")
