#!/usr/bin/env Rscript
# Optional setup; project-local library, no global package changes.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",args[grepl("^--file=",args)][1])
project <- normalizePath(file.path(dirname(script),".."))
lib <- file.path(project,"tmp","Rlib");dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
if(requireNamespace("ElPiGraph.R",quietly=TRUE) && startsWith(find.package("ElPiGraph.R"),lib)) {
  cat("Project ElPiGraph.R is already usable:",find.package("ElPiGraph.R"),"\n");quit(status=0)
}
core <- c("tictoc","irlba","ggplot2","segmented","infotheo")
for(p in core)if(!requireNamespace(p,quietly=TRUE))install.packages(p,lib=lib)
if(!requireNamespace("distutils",quietly=TRUE)) {
  if(!requireNamespace("remotes",quietly=TRUE))install.packages("remotes",lib=lib)
  remotes::install_github("Albluca/distutils",lib=lib,upgrade="never")
}
if(!requireNamespace("GenSA",quietly=TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/GenSA/GenSA_1.1.14.tar.gz",lib=lib,repos=NULL,type="source")
}
archive <- file.path(project,"tmp","ElPiGraph_official_761a48a.tar.gz")
download.file("https://api.github.com/repos/Albluca/ElPiGraph.R/tarball/761a48a",archive,mode="wb",method="libcurl")
folder <- file.path(project,"tmp","ElPiGraph_reproducible_source")
dir.create(folder,showWarnings=FALSE);untar(archive,exdir=folder)
source <- list.dirs(folder,recursive=FALSE,full.names=TRUE)[1]
description <- file.path(source,"DESCRIPTION")
file.copy(description,paste0(description,".original"),overwrite=TRUE)
d <- read.dcf(description)
imports <- trimws(strsplit(d[1,"Imports"],",")[[1]])
d[1,"Imports"] <- paste(setdiff(imports,c("ggforce","ggpubr")),collapse=", ")
write.dcf(d,description)
status <- system2(file.path(R.home("bin"),"R"),c("CMD","INSTALL",paste0("--library=",shQuote(lib)),shQuote(source)))
if(status!=0)stop("ElPiGraph.R project installation failed")
writeLines(c("Official ElPiGraph.R commit 761a48a",paste("Source:",source),
             "Only DESCRIPTION plotting Imports ggforce/ggpubr removed; R algorithm files unchanged.",
             "Plots are produced by the Python pipeline; no system font development packages required.",
             paste("GenSA compatibility pin: 1.1.14; R:",R.version.string)),file.path(lib,"StromalAtlas_setup_provenance.txt"))
cat("ElPiGraph.R installed in",lib,"\n")
