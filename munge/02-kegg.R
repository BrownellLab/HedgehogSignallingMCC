colFile <- "data/hsa04340_hedgehog_signaling_pathway.xml"

g.kegg <- parseKGML(colFile)
g <- KEGGpathway2Graph(g.kegg, genesOnly=TRUE)
outs <- sapply(edges(g), length) > 0
ins <- sapply(inEdges(g), length) > 0
ios <- outs | ins

## translate the KEGG IDs into Gene Symbol
ioGeneID <- translateKEGGID2GeneID(names(ios))
nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
names(nodesNames) <- names(ios)
    
## Get only edges between nodes
outs <- sapply(edges(hs.main), length) > 0
ins <- sapply(inEdges(hs.main), length) > 0
ios <- outs | ins

## Get a subset
keep <- nodesNames %in% hedgehog.gene.names
nodesNames <- nodesNames[keep]
ios <- ios[keep]

sg <- subGraph(names(nodesNames), g)

## To plot things, you first need to create some attributes
## Note that RCytoscape does not uses igraph format, but graph/RBGL package format
sg <- initNodeAttribute(sg, 'geneName', 'char', 'NA')
sg <- initNodeAttribute(sg, 'pathID', 'char', 'NA')
sg <- initNodeAttribute(sg, 'pathName', 'char', 'NA')
sg <- initNodeAttribute(sg, "expression", 'numeric', default.value=0)

nName <- nodes(sg)
nodeData(sg, nName, 'geneName') <- as.character(nodesNames)

new.expressionvals <- merge(data.frame(nodesNames, names=names(nodesNames)), expressionvals, by.x="nodesNames", by.y="gene.symbol")
nodeData(sg, as.character(new.expressionvals$names), 'expression') <- new.expressionvals$meanexpr
## nodeData(sg, nName, 'pathID') <- as.character(keggHS$pathID)
## nodeData(sg, nName, 'pathName') <- as.character(keggHS$pathDS)    

## ## ## you need to install and activate cytoscape plugin cytoscapeRPC to be able to send graphs directly into cytoscape

## ## Make a connection
hsCW <- new.CytoscapeWindow('sg', graph=sg)

## Create a view
displayGraph(hsCW)

setNodeLabelRule(hsCW, 'geneName')
setDefaultBackgroundColor(hsCW, '#FFFFFF')
setDefaultEdgeColor(hsCW, '#000000')
control.points <- c(0, median(expressionvals$meanexpr), max(expressionvals$meanexpr) + 1)
node.colors <- c ("#0000FF", "#0000AA", "#FFFFFF", "#FF0000", "#AA0000")
setNodeColorRule(hsCW, 'expression', control.points=control.points, colors=node.colors, mode='interpolate')
setNodeBorderColorDirect(hsCW, nName, '#000000')
setEdgeTargetArrowShapeDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.values="Arrow")
setEdgeTargetArrowColorDirect(obj=hsCW, edge.names=as.character(cy2.edge.names(hsCW@graph)), new.colors="#000000")

redraw(hsCW)

    
    
    
## sg <- subGraph(names(nodesNames), g)

## nAttrs <- list();
## nAttrs$fillcolor <- makeAttr(sg, "lightgrey", list(orange=names(ios)[ios]))

## nAttrs$label <- nodesNames
## plot(sg, "neato", nodeAttrs=nAttrs,
##      attrs=list(node=list(fillcolor="lightgreen",
##                   width="0.75", shape="ellipse"),
##        edge=list(arrowsize="0.7")))

