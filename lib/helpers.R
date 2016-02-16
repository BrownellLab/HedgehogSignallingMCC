makeAttr <- function(graph, default, valNodeList) {
  tmp <- nodes(graph)
  x <- rep(default, length(tmp)); names(x) <- tmp
  
  if(!missing(valNodeList)) {
    stopifnot(is.list(valNodeList))
    allnodes <- unlist(valNodeList)
    stopifnot(all(allnodes %in% tmp))
    for(i in seq(valNodeList)) {
      x[valNodeList[[i]]] <- names(valNodeList)[i]
    }
  }
  return(x)
}

# For KEGG to Cytoscape
process_pathway_graph <- function(pathway_graph) {
  ## Requires an item with a path_id, path_name, and a graph object.
  ## Returns a data frame with some info.
  
  pathID <- pathway_graph$path_id
  pathDes <- pathway_graph$path_name
  mapK2HS <- pathway_graph$graph

  nList <- nodes(mapK2HS)
  geneIDs <-translateKEGGID2GeneID(nList)
  namesHS <- unlist(mget(geneIDs, org.Hs.egSYMBOL,ifnotfound=NA))
  ## numHS <- unique(c(numHS, unlist(namesHS)))
  
  if (length(namesHS) > 0) {
    data.frame(HS_ID = nList, HS_gene = namesHS,
               KpathID = rep(pathID, times=length(namesHS)),
               KpathDS = rep(pathDes, times=length(namesHS)))
  }  
  
}

get_best_probe <- function(probenames) {
  ## Get the single matching probe, unless there are none
  best <- grep("[0-9]+_at", probenames)
  if (length(best) == 0) {
    probenames
  } else {
    probenames[best]
  }
}

# For KEGG to Cytoscape
## Process the data somehow....
pDS <- function(x) {data.frame(pathDS=paste(unique(x$KpathDS), collapse="|"))}
pDS <- function(x) {data.frame(pathDS=paste(unique(x$KpathDS), collapse="|"))}
pID <- function(x) {data.frame(pathID=paste(unique(x$KpathID), collapse="|"))}
# what's this for???
## pLG <- function(x) {data.frame(pathLG=paste(unique(x$KpathLG), collapse="|"))}
