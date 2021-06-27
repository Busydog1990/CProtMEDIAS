#' Calculate score of each site in alignment
#'
#' @param alignment input sequence alignment result, for DNA/RNA/protein sequences, DNAStringSet/RNAStringSet/AAStringSet or character class were required. All sequences must have the same length.
#' @param alignment_file the path of input sequence alignment result.
#' @param standard Logical. If TRUE, Only keep the score of AA_STANDARD/DNA_BASES/RNA_BASES for protein/DNA/RNA sequences. Default:FALSE.
#' @param symbol Character. Prefix of colnames of result. Default:"site".
#' @param type Input sequence type.One of "DNA", "RNA" or "AA".
#' @param seq_length_limit Integer. Lower limit of input sequence length. Default:0.
#' @param remove_duplicated Logical. If TRUE, only one sequence with repeated names in the score matrix is kept. Default:FALSE.
#' @return score matrix of each site. rownames of the repeated sequence name will be suffixed.
#' @export
#' @examples
#' require(ggseqlogo)
#' data(ggseqlogo_sample)
#' test_aa <- get_alignment_score(alignment = seqs_aa[[1]],type = "AA")
#' test_dna <- get_alignment_score(alignment = seqs_dna[[1]],type = "DNA")
get_alignment_score <- function(alignment = NULL,alignment_file = NULL,standard = F,
                                symbol = "site",type = c("DNA","RNA","AA"),
                                seq_length_limit = 0L,remove_duplicated = F){

  ##Read in the alignment sequence and convert it to matrix format
  if (!type %in% c("DNA","RNA","AA")){
    stop("type must one of 'DNA','RNA' or 'AA'")
  }

  if (is.null(alignment)){
    my_seq_alignment <- switch(type,DNA = Biostrings::readDNAStringSet(alignment_file),
                              RNA = Biostrings::readRNAStringSet(alignment_file),
                              AA = Biostrings::readAAStringSet(alignment_file))
  } else {
    if ("character" %in% class(alignment)){
      my_seq_alignment <- switch(type,DNA = Biostrings::DNAStringSet(alignment),
                                 RNA = Biostrings::RNAStringSet(alignment),
                                 AA = Biostrings::AAStringSet(alignment))
    } else if (any(class(alignment) %in% c("DNAStringSet","RNAStringSet","AAStringSet"))){
      my_seq_alignment <- alignment
    } else {
      stop("wrong alignment")
    }
  }
  stopifnot("All sequences must have the same length" =
              max(Biostrings::width(my_seq_alignment)) == min(Biostrings::width(my_seq_alignment)))
  stopifnot("All sequences must have the length above seq_length_limit" =
              max(Biostrings::width(my_seq_alignment)) >= seq_length_limit)

  if (max(Biostrings::width(my_seq_alignment)) <= 100){warning("Since length of sequences are shorter than 100,May not have enough feature sites")}

  my_alphabet <- Biostrings::alphabet(my_seq_alignment)

  my_seq_alignment_mat <- as.matrix(my_seq_alignment)

  ##Construct an amino acid percentage matrix and calculate the percentage of amino acids at each position
  my_percent <- get_seq_percent(my_seq_alignment_mat,my_alphabet)

  colnames(my_percent) <- my_alphabet
  ##Set "-", "*", "+" and "." symbol amino acid percentage to 0
  if (type == "AA"){
    my_percent[,c("-","*","+",".")] = 0
  } else {
    my_percent[,c("-","+",".")] = 0
  }
  if (standard){
    my_other <- switch(type,DNA = Biostrings::DNA_ALPHABET[!Biostrings::DNA_ALPHABET %in% Biostrings::DNA_BASES],
                       AA = Biostrings::AA_ALPHABET[!Biostrings::AA_ALPHABET %in% Biostrings::AA_STANDARD],
                       RNA = Biostrings::RNA_ALPHABET[!Biostrings::RNA_ALPHABET %in% Biostrings::RNA_BASES])
    my_percent[,my_other] = 0
  }
  ##Calculate the score of each site of each sequence according to the amino acid percentage of each site.

  my_seq_alignment_score <- matrix(0,nrow(my_seq_alignment_mat),ncol(my_seq_alignment_mat))

  rownames(my_seq_alignment_score) <- rownames(my_seq_alignment_mat)

  colnames(my_seq_alignment_score) <- paste(symbol,1:ncol(my_seq_alignment_score))

  for (i in 1:ncol(my_seq_alignment_score)){

    my_seq_alignment_score[,i] <- my_percent[i,my_seq_alignment_mat[,i]]

  }

  if (remove_duplicated){
	  my_seq_alignment_score <- my_seq_alignment_score[!duplicated(rownames(my_seq_alignment_score)),]
  }

  return(my_seq_alignment_score)
}

