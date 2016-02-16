## some basal cell carcinoma public data
## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7553
gse.basalcell <- getGEO("GSE7553", GSEMatrix = TRUE)
gse.basalcell <- gse.basalcell[['GSE7553_series_matrix.txt.gz']]

datadir.basalcell <- "/var/preserve/data/GSE7553/"
files.basalcell <- paste("/var/preserve/data/GSE7553/", sampleNames(gse.basalcell), ".CEL", sep="")

data.basalcell <- ReadAffy(filenames=files.basalcell,
                           sampleNames=sampleNames(gse.basalcell),
                           phenoData=pData(gse.basalcell))

eset.basalcell <- rma(data.basalcell)
fData(eset.basalcell) <- fData(gse.basalcell)

gse.basalcell <- gse.basalcell[, grep("Basal cell carcinoma", pData(gse.basalcell)$description)]
exprs.gse.basalcell <- log2(exprs(gse.basalcell))

## ## Get all probes by gene name, filter by best probe
## gse.basalcell.hedgehog <- subset(gse.basalcell, fData(gse.basalcell)[['Gene Symbol']] %in% hedgehog.gene.names$gene.symbol)
## best.probes <- as.character(aaply(rownames(exprs(gse.basalcell.hedgehog)), 1, get_best_probe))
## gse.basalcell.hedgehog <- gse.basalcell.hedgehog[best.probes, ]

## Get specific probes
gse.basalcell.hedgehog <- gse.basalcell[as.character(hedgehog.genes.curated$ID), ]


## some medulloblastoma public data
## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37418
gse.medullo <- getGEO("GSE37418", GSEMatrix = TRUE)
gse.medullo <- gse.medullo[['GSE37418_series_matrix.txt.gz']]

## ## Get all probes by gene name, filter by best probe
## gse.medullo.hedgehog <- subset(gse.medullo, fData(gse.medullo)[['Gene Symbol']] %in% hedgehog.gene.names$gene.symbol)
## best.probes <- as.character(aaply(rownames(exprs(gse.medullo.hedgehog)), 1, get_best_probe))
## gse.medullo.hedgehog <- gse.medullo.hedgehog[best.probes, ]

## Get specific probes
gse.medullo.hedgehog <- gse.medullo[as.character(hedgehog.genes.curated$ID), ]

## Load MCC data
load("/var/preserve/git/MCC_only/cache/eset.RData")

## ## Get all probes by gene name, filter by best probe
## eset.hedgehog <- subset(eset, fData(eset)$gene.symbol %in% hedgehog.gene.names$gene.symbol)
## best.probes <- as.character(aaply(rownames(exprs(eset.hedgehog)), 1, get_best_probe))
## eset.hedgehog <- eset.hedgehog[best.probes, ]

## Get specific probes
eset.hedgehog <- eset[as.character(hedgehog.genes.curated$ID), ]


## Common hedgehog probes
common.hedgehog.probes <- intersect(rownames(fData(eset.hedgehog)), rownames(fData(gse.basalcell.hedgehog)))

eset.hedgehog <- eset.hedgehog[common.hedgehog.probes, ]
gse.basalcell.hedgehog <- gse.basalcell.hedgehog[common.hedgehog.probes, ]
gse.medullo.hedgehog <- gse.medullo.hedgehog[common.hedgehog.probes, ]

## ## Housekeeping gene based normalization
## controls.affyid <- c("AFFX-HUMGAPDH/M33197_5_at", "AFFX-HUMGAPDH/M33197_M_at", "AFFX-HUMGAPDH/M33197_3_at",
##                      "AFFX-HSAC07/X00351_3_at", "AFFX-HSAC07/X00351_M_at", "AFFX-HSAC07/X00351_5_at")

