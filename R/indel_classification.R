
#' Classify indels into different types
#'
#' @param indels A indel list. Columns include: "Sample","chr", "position", "REF", "ALT"
#' @param genome.v : "hg19", "hg38"
#' @return Classified indel list
#' @export
indel_classifier89 <- function(indels, genome.v){
  #names(indels) <- c("Sample","chr", "position", "REF", "ALT")
  # indels$chr <- paste0("chr",indels$chr)

  # indel[indel$chr=="23","chr"]="X"
  # indel[indel$chr=="24","chr"]="Y"

  ## prepare indels
  prep_df <- prepare_indels(as.data.frame(indels),"pancan",genome.v)
  ## Segment indels
  s <- segment_indels(prep_df) #get('segment_indels',envir = .GlobalEnv)(prep_df$change,prep_df$slice3)
  s$mh_length <- s$indel.length-s$spacer_length

  ## Assign indels m5
  s_classified <- assign_channels_m5(s)

  return(s_classified)

}
