% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_monocle_workflow.R
\name{site_score}
\alias{site_score}
\title{Calculate the match score between two sites.}
\usage{
site_score(site1, site2, substitutionMatrix, quickly = F)
}
\arguments{
\item{site1}{score of site1.}

\item{site2}{score of site2.}

\item{substitutionMatrix}{Substitution matrix representing the fixed substitution scores for an alignment.}

\item{quickly}{If TRUE, only keep the score of AA_STANDARD/DNA_BASES/RNA_BASES for protein/DNA/RNA sequences. Default:TRUE.}
}
\value{
Combined multiple alignment results.
}
\description{
Calculate the match score between two sites.
}
