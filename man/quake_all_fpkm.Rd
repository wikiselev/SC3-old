% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/data-examples.R
\docType{data}
\name{quake_all_fpkm}
\alias{quake_all_fpkm}
\title{Single cell RNA-Seq data extracted from a publication by Steven Quake.}
\format{\preformatted{ num [1:23271, 1:80] 0 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:23271] "0610005C13Rik" "0610007C21Rik" "0610007L01Rik" "0610007N19Rik" ...
  ..$ : chr [1:80] "2" "4" "3" "2" ...
}}
\usage{
quake_all_fpkm
}
\description{
Source: Treutlein, B. et al. Reconstructing lineage hierarchies of the distal
lung epithelium using single-cell RNA-seq. Nature 509, 371–375 (2014).
}
\details{
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52583}

Columns represent cells, rows represent genes expression values. Colnames
respresent indexes of cell clusters (known information based on the
experimental protocol). There are 80 cells and 5 clusters in this dataset.
}
\keyword{datasets}

