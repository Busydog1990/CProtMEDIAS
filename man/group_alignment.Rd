% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_monocle_workflow.R
\name{group_alignment}
\alias{group_alignment}
\title{Combine two alignment results.}
\usage{
group_alignment(
  alignment1,
  alignment2,
  seqtype = "AA",
  gapOpening = -4,
  gapExtension = -1,
  substitutionMatrix = "BLOSUM62",
  standard = T,
  alignment_thres = 0
)
}
\arguments{
\item{alignment1}{Alignment result of group1.}

\item{alignment2}{Alignment result of group2.}

\item{seqtype}{One of DNA, RNA and AA.}

\item{gapOpening}{Numeric. The cost for opening a gap in the alignment.}

\item{gapExtension}{Numeric. The incremental cost incurred along the length of the gap in the alignment.}

\item{substitutionMatrix}{Substitution matrix representing the fixed substitution scores for an alignment.}

\item{standard}{Logical. If TRUE, only keep the score of AA_STANDARD/DNA_BASES/RNA_BASES for protein/DNA/RNA sequences. Default:TRUE.}

\item{alignment_thres}{Numeric. Remove sites with amino acid ratio lower than alignment_thres. Default:0.}
}
\value{
Combined multiple alignment results.
}
\description{
Solves (Needleman-Wunsch) global alignment problems on two alignment results.
}
