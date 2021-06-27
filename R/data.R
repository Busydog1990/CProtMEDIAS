#' Multiple sequences alignment of homeobox transcription factor superfamily
#'
#' A dataset containing first 200 sequences in multiple sequences
#' alignment result of DNA binding domain (DBD) of HB_other, HB_PHD,
#' HD_ZIP, TALE and WOX transcription factor family
#' @format A list with 5 elements:
#' \describe{
#'   \item{HB_other}{multiple sequences alignment result of DNA
#'                   binding domain (DBD) of HB_other,
#'      \url{http://planttfdb.gao-lab.org/multi_align/family/HB-other/domain_aln.fas}}
#'   \item{HB_PHD}{multiple sequences alignment result of DNA
#'                   binding domain (DBD) of HB_PHD,
#'      \url{http://planttfdb.gao-lab.org/multi_align/family/HB-PHD/domain_aln.fas}}
#'   \item{HD_ZIP}{multiple sequences alignment result of DNA
#'                   binding domain (DBD) of HD_ZIP,
#'      \url{http://planttfdb.gao-lab.org/multi_align/family/HD_ZIP/domain_aln.fas}}
#'   \item{TALE}{multiple sequences alignment result of DNA
#'                   binding domain (DBD) of TALE,
#'      \url{http://planttfdb.gao-lab.org/multi_align/family/TALE/domain_aln.fas}}
#'   \item{WOX}{multiple sequences alignment result of DNA
#'                   binding domain (DBD) of WOX,
#'      \url{http://planttfdb.gao-lab.org/multi_align/family/WOX/domain_aln.fas}}
#' }
"Homeobox_align"


#' Annotation of DNA binding domain of homeobox transcription factor superfamily
#'
#'
#' @format A data frame with 1000 rows and 11 variables:
#' \describe{
#'   \item{Species}{Species of transcription factors}
#'   \item{ID}{Gene id of transcription factors}
#'   \item{Names}{Names of DNA binding domain}
#'   \item{Family}{Family of transcription factors}
#'   \item{Description}{Description of transcription factor family}
#'   \item{Taxonomic1}{Phylum corresponding to each species}
#'   \item{Taxonomic2}{Phylum or Subclass corresponding to each species}
#'   \item{Alias}{Alias of species}
#'   \item{Genome}{Species with or without genome sequence,
#'     \url{http://planttfdb.gao-lab.org/index.php}}
#'   \item{WOX_subgroup}{Subfamily of WOX transcription factors family}
#'   \item{At_WOX}{The most similar WOX gene in Arabidopsis, only for TF of WOX family}
#' }
#'
"Homeobox_annot"

