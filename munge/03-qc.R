## Do quality control on batch corrected data
qc.all.batcheffect <- arrayQualityMetrics(expressionset=neweset.batcheffect, outdir="reports/QC/RMA.batcheffect.all/", force=T, intgroup="cancertype")
