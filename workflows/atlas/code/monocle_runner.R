#!/usr/bin/env Rscript
# Monocle 2 DDRTree on mean log-expression meta-cells; never invent raw counts.
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]; variant <- args[2]; output <- args[3]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(jsonlite))
set.seed(42)

# Process-local compatibility with igraph >=2; installed libraries are untouched.
ns <- asNamespace("monocle")
for (nm in ls(ns, all.names=TRUE)) {
  obj <- get(nm, envir=ns)
  if (!is.function(obj)) next
  txt <- paste(deparse(obj), collapse="\n")
  revised <- gsub("\\bnei\\(", ".nei(", txt, perl=TRUE)
  revised <- gsub("neimode =", "mode =", revised, fixed=TRUE)
  revised <- gsub('if (class(projection) != "matrix")', 'if (!is.matrix(projection))', revised, fixed=TRUE)
  if (!identical(txt,revised)) {
    replacement <- eval(parse(text=revised), envir=environment(obj))
    unlockBinding(nm,ns);assign(nm,replacement,envir=ns);lockBinding(nm,ns)
  }
}
spec <- fromJSON(file.path(input,"prepared.json"),simplifyVector=FALSE)
obs <- read.csv(file.path(input,"metadata.csv"),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
genes <- read.csv(file.path(input,"genes.csv"),stringsAsFactors=FALSE)
con <- file(file.path(input,"expression.f32"),"rb")
vals <- readBin(con, what="numeric", n=as.double(spec$n_cells)*spec$n_genes, size=4, endian="little")
close(con)
expr <- matrix(vals,nrow=spec$n_genes,ncol=spec$n_cells)
rm(vals);gc()
rownames(expr) <- genes$gene;colnames(expr) <- rownames(obs)
fd <- data.frame(gene_short_name=genes$gene,row.names=genes$gene)
cds <- newCellDataSet(expr,phenoData=new("AnnotatedDataFrame",data=obs),
                     featureData=new("AnnotatedDataFrame",data=fd),
                     expressionFamily=VGAM::uninormal(),lowerDetectionLimit=0)
pData(cds)$Size_Factor <- 1
ordering <- unlist(spec$variants[[variant]],use.names=FALSE)
cds <- setOrderingFilter(cds,ordering)
cat("DDRTree_START",variant,ncol(cds),"meta-cells",format(Sys.time()),"\n");flush.console()
cds <- reduceDimension(cds,max_components=2,reduction_method="DDRTree",norm_method="none",
                       relative_expr=FALSE,pseudo_expr=0,scaling=TRUE,ncenter=100,maxIter=30,
                       verbose=FALSE)
cat("DDRTree_DONE",format(Sys.time()),"\n");flush.console()
principal_nodes <- reducedDimK(cds)
principal_tree <- minSpanningTree(cds)
cds <- orderCells(cds)
cat("ORDERING_DONE",format(Sys.time()),"\n");flush.console()
state <- as.character(pData(cds)$State)
root_type <- spec$root_state
tab <- table(state, pData(cds)$general==root_type)
root_n <- if("TRUE" %in% colnames(tab)) tab[,"TRUE"] else rep(0,nrow(tab))
enrich <- root_n/rowSums(tab)
enrich[root_n < max(5,0.01*sum(root_n))] <- -Inf
root_state <- rownames(tab)[which.max(enrich)]
cds <- orderCells(cds,root_state=root_state)
xy <- t(reducedDimS(cds))
result <- data.frame(cell_id=rownames(obs),x=xy[,1],y=xy[,2],pseudotime=pData(cds)$Pseudotime,
                     branch=as.character(pData(cds)$State),general=obs$general)
write.csv(result,file.path(output,"cells.csv"),row.names=FALSE)
nodes <- data.frame(node=colnames(principal_nodes),x=principal_nodes[1,],y=principal_nodes[2,])
if(is.null(nodes$node))nodes$node <- as.character(seq_len(ncol(principal_nodes)))
write.csv(nodes,file.path(output,"nodes.csv"),row.names=FALSE)
ed <- as.data.frame(igraph::as_edgelist(principal_tree,names=TRUE),stringsAsFactors=FALSE)
colnames(ed) <- c("source","target");write.csv(ed,file.path(output,"edges.csv"),row.names=FALSE)
saveRDS(cds,file.path(output,"monocle_cds.rds"))
write_json(list(method="Monocle2-DDRTree",version=as.character(packageVersion("monocle")),
                ordering_variant=variant,n_genes=length(ordering),root_annotation=root_type,
                root_monocle_state=root_state,expression="mean log-expression; uninormal family, no re-log, no count likelihood",
                compatibility="Process-local nei/neimode/project2MST adaptations for igraph >=2"),
           file.path(output,"method.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
cat("MONOCLE_SUCCESS",variant,nrow(result),"meta-cells\n")
