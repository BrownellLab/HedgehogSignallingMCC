hedgehog.gene.names <- data.frame(gene.symbol=c("SHH", "WNT1", "SMO", "PTCH1", "SUFU", "HHIP", "GAS1", "ZIC2", "BMP4", "GLI1", "GLI2", "GLI3", "BTRC", "LRP2"))

## These were picked by looking at the medulloblastoma data specifically for the SHH subgroup.
## The PTCH1 gene has one probe that is highly expressed in all samples (even MCC, non-SHH medullos)
hedgehog.genes.curated <- data.frame(ID=c("208570_at", "218629_at", "207586_at", "209816_at", "205710_at", "205201_at", "208057_s_at",
                                       "207034_s_at", "206646_at", "204457_s_at", "204456_s_at", "222374_at", "204901_at", "211518_s_at"),
                                     gene.symbol=c("WNT1", "SMO", "SHH", "PTCH1", "LRP2", "GLI3", "GLI2", "GLI2", "GLI1", "GAS1", "GAS1", "BTRC", "BTRC", "BMP4"))
