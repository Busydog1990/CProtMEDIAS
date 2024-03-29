% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_monocle_workflow.R
\name{get_seq_conservation}
\alias{get_seq_conservation}
\title{Get sequence conservation of each site.}
\usage{
get_seq_conservation(
  seqs = NULL,
  alignment_file = NULL,
  type = c("DNA", "RNA", "AA"),
  group = NULL,
  ...
)
}
\arguments{
\item{seqs}{input sequence alignment result, for DNA/RNA/protein sequences, DNAStringSet/RNAStringSet/AAStringSet or character class were required. All sequences must have the same length.}

\item{alignment_file}{The path of input sequence alignment result.}

\item{type}{Character. One of DNA, RNA and AA.}

\item{group}{Vector. Grouping information of each sequence. The length of group vector is consistent with the number of sequences.}

\item{...}{Other parameters of \code{\link{consensusMatrix}}.}
}
\value{
sequence conservation matrix of each site.
}
\description{
Get sequence conservation of each site.
}
