expr.data <- data.frame(exprs(gse.basalcell.hedgehog), ID=rownames(exprs(gse.basalcell.hedgehog)), gene.symbol=fData(gse.basalcell.hedgehog)[['Gene Symbol']])

## Now set expression values
expressionvals <- ddply(expr.data, "gene.symbol", function(x) apply(x[, colnames(exprs(gse.basalcell.hedgehog))], 2, mean))
expressionvals <- transform(expressionvals, meanexpr=apply(expressionvals[, colnames(exprs(gse.basalcell.hedgehog))], 1, mean))
new.expressionvals <- merge(data.frame(nodesNames, names=names(nodesNames)), expressionvals, by.x="nodesNames", by.y="gene.symbol")

sg.basal <- subGraph(names(nodesNames), g)
sg.basal <- initNodeAttribute(sg.basal, 'geneName', 'char', 'NA')
sg.basal <- initNodeAttribute(sg.basal, "expression", 'numeric', default.value=0)

nodeData(sg.basal, as.character(new.expressionvals$names), 'expression') <- new.expressionvals$meanexpr

nName <- nodes(sg.basal)
nodeData(sg.basal, nName, 'geneName') <- as.character(nodesNames)

## Create a view
hsCW <- new.CytoscapeWindow('sg.basal', graph=sg.basal)
displayGraph(hsCW)

setNodeLabelRule(hsCW, 'geneName')
setDefaultBackgroundColor(hsCW, '#FFFFFF')
setDefaultEdgeColor(hsCW, '#000000')
control.points <- c(0, mean(exprs(gse.basalcell)), max(exprs(gse.basalcell)))
node.colors <- c ("#0000FF", "#0000FF", "#FFFFFF", "#FF0000", "#FF0000")
setNodeColorRule(hsCW, 'expression', control.points=control.points, colors=node.colors, mode='interpolate')
setNodeBorderColorDirect(hsCW, nName, '#000000')
setEdgeTargetArrowShapeDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.values="Arrow")
setEdgeTargetArrowColorDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.colors="#000000")

redraw(hsCW)

