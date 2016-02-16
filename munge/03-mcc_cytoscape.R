expr.data <- data.frame(exprs(eset.hedgehog), ID=rownames(exprs(eset.hedgehog)), gene.symbol=fData(eset.hedgehog)$gene.symbol)

sg.mcc <- subGraph(names(nodesNames), g)
sg.mcc <- initNodeAttribute(sg.mcc, 'geneName', 'char', 'NA')
sg.mcc <- initNodeAttribute(sg.mcc, "expression", 'numeric', default.value=0)

## Now set expression values

expressionvals <- ddply(expr.data, "gene.symbol", function(x) apply(x[, colnames(exprs(eset.hedgehog))], 2, mean))
expressionvals <- transform(expressionvals, meanexpr=apply(expressionvals[, colnames(exprs(eset.hedgehog))], 1, mean))
new.expressionvals <- merge(data.frame(nodesNames, names=names(nodesNames)), expressionvals, by.x="nodesNames", by.y="gene.symbol")
nodeData(sg.mcc, as.character(new.expressionvals$names), 'expression') <- new.expressionvals$meanexpr

nName <- nodes(sg.mcc)
nodeData(sg.mcc, nName, 'geneName') <- as.character(nodesNames)

## Create a view
hsCW <- new.CytoscapeWindow('sg.mcc', graph=sg.mcc)
displayGraph(hsCW)

setNodeLabelRule(hsCW, 'geneName')
setDefaultBackgroundColor(hsCW, '#FFFFFF')
setDefaultEdgeColor(hsCW, '#000000')
control.points <- c(0, mean(exprs(eset)), max(exprs(eset)))
node.colors <- c ("#0000FF", "#0000FF", "#FFFFFF", "#FF0000", "#FF0000")
setNodeColorRule(hsCW, 'expression', control.points=control.points, colors=node.colors, mode='interpolate')
setNodeBorderColorDirect(hsCW, nName, '#000000')
setEdgeTargetArrowShapeDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.values="Arrow")
setEdgeTargetArrowColorDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.colors="#000000")

redraw(hsCW)

    
    
    
## sg.mcc <- subGraph(names(nodesNames), g)

## nAttrs <- list();
## nAttrs$fillcolor <- makeAttr(sg.mcc, "lightgrey", list(orange=names(ios)[ios]))

## nAttrs$label <- nodesNames
## plot(sg.mcc, "neato", nodeAttrs=nAttrs,
##      attrs=list(node=list(fillcolor="lightgreen",
##                   width="0.75", shape="ellipse"),
##        edge=list(arrowsize="0.7")))

