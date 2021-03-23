# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# Plot indel profile/signature in a 89-channel bar plot for single sample
gen_plot_catalouge89_single<- function(muts_basis,text_size,plot_title){
  #mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  indel_template_type_4 <- data.frame("IndelType"=c("A|[+C]Rep=0|A","A|[+C]Rep=0|T","[+C]Rep_leq3","[+C]Rep_456","[+C]Rep_789",


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
                                                    "Complex"),

                                      "Indel"=c("[+C]","[+C]","[+C]","[+C]","[+C]",


                                                "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",
                                                "[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]","[+T]",

                                                "Ins_NonRep","Ins_NonRep","Ins_NonRep","Ins_NonRep",
                                                "Ins_nMer","Ins_nMer",

                                                "[-C]","[-C]","[-C]","[-C]","[-C]","[-C]","[-C]","[-C]",
                                                "[-C]","[-C]",


                                                "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",
                                                "[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]","[-T]",


                                                "Del_NonRep","Del_NonRep",
                                                "Del_nMer","Del_nMer","Del_nMer","Del_nMer",
                                                "Del_Spaced","Del_Spaced","Del_Spaced","Del_Spaced","Del_Spaced","Del_Spaced","Del_Spaced",
                                                "Complex")
  )

  muts_basis_melt <- melt(muts_basis,"IndelType")

  muts_basis_melt <- merge(indel_template_type_4, muts_basis_melt,by="IndelType",all.x=T)
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c("IndelType","Indel","Sample","freq")
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)


  # indel_mypalette_fill <- c("skyblue","orange", "blue","tomato","greenyellow","pink","grey","purple","deeppink","black")
  indel_mypalette_fill <- c("#06C2F4", # [-C]
                            "#FF8642", # [-T]
                            "#007996", # [+C]
                            "#C0362C", # [+T]
                            "#FEABB9", # Complex
                            "#B1DDA1", # Del_nMer
                            "#C3B7AC", # Del_NonRep
                            "#91278F", # Del_Spaced
                            "#668D3C", # Ins_nMer
                            "#816C5B")   # Ins_NonRep

  indel_positions <- indel_template_type_4$IndelType
  entry <- table(indel_template_type_4$Indel)
  order_entry <- c("[+C]", "[+T]", "Ins_NonRep", "Ins_nMer", "[-C]", "[-T]", "Del_NonRep", "Del_nMer", "Del_Spaced", "Complex")
  entry <- entry[order_entry]
  blocks <- data.frame(Type=unique(indel_template_type_4$Indel),
                       fill=indel_mypalette_fill,
                       xmin=c(0,cumsum(entry)[-length(entry)])+0.5,
                       xmax=cumsum(entry)+0.5)
  blocks$ymin <- max(muts_basis_melt$freq)*1.08
  blocks$ymax <- max(muts_basis_melt$freq)*1.2
  blocks$labels <-c("+C", "+T", "+N", "+M", "-C", "-T", "-N", "-M", "-Mh", "X")


  p <- ggplot(data=muts_basis_melt, aes(x=IndelType, y=freq,fill=Indel))+ geom_bar(stat="identity",position="dodge", width=.7)+xlab("Indel Types")+ylab("Count")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions)+ggtitle(plot_title)
  p <- p+scale_fill_manual(values=indel_mypalette_fill)
  p <- p+theme_classic()+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=5,colour = "black",hjust=1),
                               axis.text.y=element_text(size=10,colour = "black"),
                               axis.line.y=element_blank(),
                               legend.position = "none",
                               axis.title.x = element_text(size=15),
                               axis.title.y = element_text(size=15))
  ## Add the overhead blocks

  p <- p+geom_rect(data = blocks, aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=Type),inherit.aes = F)+
    # geom_text(data=blocks,aes(x=(xmax+xmin)/2,y=(ymax+ymin)/2,label=labels),size=text_size,inherit.aes = F,colour="white")
    geom_text(data=blocks,aes(x=(xmax+xmin)/2,y=ymax*1.1,label=labels),size=text_size,inherit.aes = F)

  return(p)


  # return(muts_basis)
}


# Plot indel profile/signature in a 89-channel bar plot, original plots_type_4_m4_89 function
plots_indel_89ch<- function(muts_basis,colnum, h,w,text_size,outputname){

  muts_basis2 <- muts_basis[,names(muts_basis) != "IndelType"]
  p_all <- list()
  for(i in 1:dim(muts_basis2)[2]){

    p <- gen_plot_catalouge89_single(data.frame("Sample"=muts_basis2[,i],"IndelType"=rownames(muts_basis2)), text_size,names(muts_basis2)[i])
    p_all[[length(p_all)+1]] <- p

  }

  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)

  do.call("grid.arrange", c(p_all, ncol = colnum))


  dev.off()

}
