#' Generate indel catalogue in 89 channels
#'
#' @param muts_list A indel list
#' @param sample_col Sample column name
#' @return A 89 channel indel catalogue
#' @export
gen_catalogue89 <- function(muts_list, sample_col){
  indel_catalogue <- data.frame(table(muts_list[,sample_col],muts_list$type_4))
  names(indel_catalogue) <- c("Sample","type_4","freq")
  indel_catalogue <- reshape2::dcast(indel_catalogue,type_4~Sample,value.var="freq")

  indel_template_type_4 <- data.frame("type_4"=c("A|[+C]Rep=0|A","A|[+C]Rep=0|T","[+C]Rep_leq3","[+C]Rep_456","[+C]Rep_789",


                                                 "A|[+T]Rep_leq4|A","A|[+T]Rep_leq4|C","A|[+T]Rep_leq4|G","C|[+T]Rep_leq4|A","C|[+T]Rep_leq4|C","C|[+T]Rep_leq4|G","G|[+T]Rep_leq4|A","G|[+T]Rep_leq4|C","G|[+T]Rep_leq4|G",
                                                 "A|[+T]Rep_567|A","A|[+T]Rep_567|C","A|[+T]Rep_567|G","C|[+T]Rep_567|A","C|[+T]Rep_567|C","C|[+T]Rep_567|G","G|[+T]Rep_567|A","G|[+T]Rep_567|C","G|[+T]Rep_567|G",
                                                 "A|[+T]Rep_89|A","A|[+T]Rep_89|C","A|[+T]Rep_89|G","C|[+T]Rep_89|A","C|[+T]Rep_89|C","C|[+T]Rep_89|G","G|[+T]Rep_89|A","G|[+T]Rep_89|C","G|[+T]Rep_89|G",

                                                 "Ins_NonRep_R0_L234","Ins_NonRep_R0_L5","Ins_NonRep_R1_L234","Ins_NonRep_R1_L5",

                                                 "Ins_nMer_R234","Ins_nMer_R5",

                                                 "[-C]Rep=1|A","[-C]Rep=1|T",
                                                 "[-C]Rep=2|A","[-C]Rep=2|T",
                                                 "[-C]Rep=3|A","[-C]Rep=3|T",
                                                 "[-C]Rep_45|A","[-C]Rep_45|T",
                                                 "[-C]Rep_leq5|G",
                                                 "[-C]Rep_6",

                                                 "A|[-T]Rep_leq4|A","A|[-T]Rep_leq4|C","A|[-T]Rep_leq4|G","C|[-T]Rep_leq4|A","C|[-T]Rep_leq4|C","C|[-T]Rep_leq4|G","G|[-T]Rep_leq4|A","G|[-T]Rep_leq4|C","G|[-T]Rep_leq4|G",
                                                 "A|[-T]Rep_567|A","A|[-T]Rep_567|C","A|[-T]Rep_567|G","C|[-T]Rep_567|A","C|[-T]Rep_567|C","C|[-T]Rep_567|G","G|[-T]Rep_567|A","G|[-T]Rep_567|C","G|[-T]Rep_567|G",
                                                 "A|[-T]Rep_89|A","A|[-T]Rep_89|C","A|[-T]Rep_89|G","C|[-T]Rep_89|A","C|[-T]Rep_89|C","C|[-T]Rep_89|G","G|[-T]Rep_89|A","G|[-T]Rep_89|C","G|[-T]Rep_89|G",

                                                 "Del_NonRep_L234","Del_NonRep_L5",
                                                 "Del_nMer_U12_R234","Del_nMer_U12_R5","Del_nMer_U3_R2","Del_nMer_U3_R3",
                                                 "Del_Spaced_short_leq5_mh1","Del_Spaced_short_leq5_mh2","Del_Spaced_short_leq5_mh3",
                                                 "Del_Spaced_long_g5_mh1","Del_Spaced_long_g5_mh2", "Del_Spaced_long_g5_mh3","Del_Spaced_long_g5_mh4",
                                                 "Complex"))


  indel_catalogue <- merge(indel_template_type_4,indel_catalogue,by="type_4",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  rownames(indel_catalogue) <- indel_catalogue[,"type_4"]
  return(indel_catalogue[,-1])
}


#' Generate indel catalogue in 21 channels
#'
#' @param muts_list A indel list
#' @param sample_col Sample column name
#' @return A 21 channel indel catalogue
#' @export
gen_catalogue21 <- function(muts_list, sample_col){
  indel_catalogue <- data.frame(table(muts_list[,sample_col],muts_list$type_3))
  names(indel_catalogue) <- c("Sample","type_3","freq")
  indel_catalogue <- reshape2::dcast(indel_catalogue,type_3~Sample,value.var="freq")

  indel_template_type_3 <- data.frame("type_3"=c("[+C]NonRep","[+C]ShortRep_leq4","[+C]LongRep_g4","[+T]NonRep","[+T]ShortRep_leq4","[+T]LongRep_g4","Ins_NonRep","Ins_nMer_ShortRep_leq4","Ins_nMer_LongRep_g4",
                                                 "[-C]NonRep","[-C]ShortRep_leq4","[-C]LongRep_g4","[-T]NonRep","[-T]ShortRep_leq4","[-T]LongRep_g4","Del_NonRep","Del_nMer_ShortRep_leq4","Del_nMer_LongRep_g4","Del_Spaced_short_leq5","Del_Spaced_long_g5",
                                                 "Complex"))
  indel_catalogue <- merge(indel_template_type_3,indel_catalogue,by="type_3",all.x=T)
  indel_catalogue[is.na(indel_catalogue)] <- 0
  rownames(indel_catalogue) <- indel_catalogue[,"type_3"]
  return(indel_catalogue[,-1])

}




