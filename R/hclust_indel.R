#' Separate indel catlogue according to a indel number cutoff
#'
#' @param mut_cat An indel catalogue of multiple samples
#' @param cutoff Threshold of indel number to separate samples according to indel burden
#' @param outputname file names
#' @return An list of two data.frame
#' @export
divid_indelcat <- function(mut_cat,cutoff=500,outputname){
  ogi_burden <- data.frame("Sample"=colnames(mut_cat),"Freq"=colSums(mut_cat))
  ogi_burden <- ogi_burden[order(ogi_burden$Freq, decreasing = T),]

  # mut_cat_lb: indel catlogue of samples with low mutation burden (< cutoff)
  mut_cat_lb <- mut_cat[,colnames(mut_cat)%in% ogi_burden[ogi_burden$Freq<cutoff,"Sample"]]

  # mut_cat_hb: indel catlogue of samples with high mutation burden (>= cutoff)
  mut_cat_hb <- mut_cat[,colnames(mut_cat)%in% ogi_burden[ogi_burden$Freq>=cutoff,"Sample"]]

  write.table(mut_cat_lb, paste0(outputname,"_cat_less",cutoff,".txt"), sep = "\t", col.names = T, row.names = T, quote = F)
  write.table(mut_cat_hb, paste0(outputname,"_cat_geq",cutoff,".txt"), sep = "\t", col.names = T, row.names = T, quote = F)

  divid_cat <- list()
  divid_cat$lb_cat <- mut_cat_lb
  divid_cat$hb_cat <- mut_cat_hb

  return(divid_cat)

}



#' Hierarchical Clustering to identify high-burden clusters
#'
#' @param mut_cat An indel catalogue of multiple samples
#' @param hclust_method This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param tissue_type Tissue type
#' @param ncluster An integer scalar or vector with the desired number of clusters
#' @param cutoff_height Numeric scalar or vector with heights where the tree should be cut
#' @return An object of class hclust which describes the tree produced by the clustering process
#' @export
hclust_indel<- function(mut_cat,hclust_method = "complete",tissue_type,ncluster=6, cutoff_height=NULL){

  indel_template_type_4_m5 <- data.frame("type_4"=c("A|[+C]Rep=0|A","A|[+C]Rep=0|T","[+C]Rep_leq3","[+C]Rep_456","[+C]Rep_789",


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


  mut_cat <- mut_cat[match(indel_template_type_4_m5$type_4, rownames(mut_cat)),]

  mut_cat <- t(mut_cat)

  # wth <- 10+round(dim(mut_cat_normal)[1]/40)
  #  ht <- w/2
  # Hierarchical clustering

  my_palette <- grDevices::colorRampPalette(c("white", "red"))(n = 100)

  # Plot the obtained dendrogram
  pdf(file=paste0(tissue_type,"_", hclust_method,"_heatmap2_normal.pdf"), onefile=TRUE,width=20,height=20, useDingbats=FALSE)
  gplots::heatmap.2(mut_cat,
                    distfun = function(x) stats::dist(x, method="euclidean"),
                    hclustfun=function(x) stats::hclust(x,method = hclust_method),
                    reorderfun=function(d, w) stats::reorder(d, w, agglo.FUN = mean), # Reorder dendrogram by branch means rather than sums
                    main = "Heatmap of indel profiles", # heat map title
                    notecol="black",      # change font color of cell labels to black
                    density.info="none",  # turns off density plot inside color legend
                    trace="none",         # turns off trace lines inside the heat map
                    margins =c(12,9),     # widens margins around plot
                    keysize=1,
                    col=my_palette,       # use on color palette defined earlier
                    key = TRUE,
                    #      breaks=col_breaks,    # enable color transition at specified limits
                    dendrogram="row",     # only draw a row dendrogram
                    Colv="NA")            # turn off column clustering
  dev.off()


  d <- stats::dist(mut_cat, method = "euclidean")
  # cossim <- 1-d^2/2
  # Hierarchical clustering using Complete Linkage
  hc1 <- stats::hclust(d, method = hclust_method )
  if(is.null(cutoff_height)==FALSE & is.null(ncluster)==TRUE){
    sub_grp <- dendextend::cutree(hc1, h = cutoff_height)
    ncluster <- length(table(sub_grp))
  }

  # avg_dend_obj <- as.dendrogram(hc1)
  # avg_col_dend <- color_branches(avg_dend_obj, h = cutoff_height)

  pdf(file=paste0(tissue_type,"_", hclust_method,"_hcluster_color.pdf"), onefile=TRUE,width=20,height=20, useDingbats=FALSE)
  # Plot the obtained dendrogram
  plot(hc1, cex = 0.6, hang = -1)
  stats::rect.hclust(hc1, k = ncluster)
  #  plot(avg_col_dend)
  dev.off()



  return(hc1)

}
