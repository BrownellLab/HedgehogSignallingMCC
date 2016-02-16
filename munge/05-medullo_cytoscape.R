## take only those that are SHH driven subgroup
gse.medullo.hedgehog.SHH <- gse.medullo.hedgehog[, pData(gse.medullo.hedgehog)$characteristics_ch1 == "subgroup: SHH"]

expr.data <- data.frame(exprs(gse.medullo.hedgehog.SHH), ID=rownames(exprs(gse.medullo.hedgehog.SHH)), gene.symbol=fData(gse.medullo.hedgehog.SHH)[['Gene Symbol']])

## Now set expression values
expressionvals <- ddply(expr.data, "gene.symbol", function(x) apply(x[, colnames(exprs(gse.medullo.hedgehog.SHH))], 2, mean))
expressionvals <- transform(expressionvals, meanexpr=apply(expressionvals[, colnames(exprs(gse.medullo.hedgehog.SHH))], 1, mean))
new.expressionvals <- merge(data.frame(nodesNames, names=names(nodesNames)), expressionvals, by.x="nodesNames", by.y="gene.symbol")

sg.medullo <- subGraph(names(nodesNames), g)
sg.medullo <- initNodeAttribute(sg.medullo, 'geneName', 'char', 'NA')
sg.medullo <- initNodeAttribute(sg.medullo, "expression", 'numeric', default.value=0)

nodeData(sg.medullo, as.character(new.expressionvals$names), 'expression') <- new.expressionvals$meanexpr

nName <- nodes(sg.medullo)
nodeData(sg.medullo, nName, 'geneName') <- as.character(nodesNames)

## Create a view
hsCW <- new.CytoscapeWindow('sg.medullo', graph=sg.medullo)
displayGraph(hsCW)

setNodeLabelRule(hsCW, 'geneName')
setDefaultBackgroundColor(hsCW, '#FFFFFF')
setDefaultEdgeColor(hsCW, '#000000')
control.points <- c(0, mean(exprs(gse.medullo)), max(exprs(gse.medullo)))
node.colors <- c ("#0000FF", "#0000FF", "#FFFFFF", "#FF0000", "#FF0000")
setNodeColorRule(hsCW, 'expression', control.points=control.points, colors=node.colors, mode='interpolate')
setNodeBorderColorDirect(hsCW, nName, '#000000')
setEdgeTargetArrowShapeDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.values="Arrow")
setEdgeTargetArrowColorDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.colors="#000000")

redraw(hsCW)

