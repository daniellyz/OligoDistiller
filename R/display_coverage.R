#' Display sequence coverage
#'
#' The function returns a ggplot of labelled fragments and internal fragments on an oligonucleotide sequence
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous geom_segment aes geom_segment geom_text labs theme_minimal theme element_blank
#' @importFrom stringr str_remove_all
#' @export

display_coverage<-function(scan.deconvoluted.annotated, seq, int.frag = F, return.plot = T){
  
# seq = "OH-Gd-Ad-Gd-Ad-Td-Cd-Td-Cd-Td-Gd-Cd-Td-Td-Cd-Td-Gd-Ad-Td-Gd-Gd-Cd-Td-Cd-Td-Cd-Td-Gd-Gd-Td-Td-Ad-Cd-Td-Gd-Cd-Cd-Ad-Gd-Td-Td-Gd-Ad-Ad-Td-Cd-Td-Gd-OH"
  
  options(stringsAsFactors = F)
  gpl = gpk = NULL
  rc1 = 0
  rc2 = 0
  bv = 0
  
  tmp_elements = strsplit(seq, "-")[[1]]
  tmp_elements = tmp_elements[2:(length(tmp_elements)-1)]
  tmp_elements = sapply(tmp_elements, function(x) str_remove_all(x, "[fdm\\*]"))
  tmp_elements = as.character(tmp_elements)
  NL = length(tmp_elements)
  seq1 = paste0(tmp_elements, collapse  = "")

  features = scan.deconvoluted.annotated$feature
  features = features[features$FORMULA!="Unknown",]
  features = features[features$FEATURE!="FLP",]
  features = features[order(features$FEATURE),]
  
  features1 = features[which(!str_detect(features$FEATURE, "-")),,drop = F]
  fragments  = features1$FEATURE
  features2 = features[which(str_detect(features$FEATURE, "-")),,drop = F]
  internals  = features2$FEATURE
  
  if (NL>2 & nrow(features1)>0){
  
    text_left = text_right = rep("", NL)
    cleavage_sites = rep("Not Cleavage Site", NL)
    cleavage_sites1 = rep("Not Cleavage Site", NL)
    cleavage_sites2 = rep("Not Cleavage Site", NL)

    for (i in 1:length(fragments)){
  
      cp = gsub("[[:digit:]]","",fragments[i])
      ft = gsub("[^0-9.-]", "", fragments[i])
      ft = as.numeric(ft)
      if (cp %in% c("a", "b", "c", "d", "aB", "bB", "cB", "dB")){
        text_left[ft] = paste0(c(text_left[ft], fragments[i]), collapse = "\n")
        cleavage_sites[ft] = "Cleavage Site"
        cleavage_sites1[ft] = "Cleavage Site left"
      } else {
        ft1 = NL - ft
        text_right[ft1] = paste0(c(text_right[ft1], fragments[i]), collapse = "\n")
        cleavage_sites[ft1] = "Cleavage Site"
        cleavage_sites2[ft1] = "Cleavage Site right"
      }
    }

  # recovery rate
    
    tmp_text = sapply(paste0(text_left, text_right), function(x) length(strsplit(x, "\n")[[1]])-1)
    rc1 = round(sum(cleavage_sites == "Cleavage Site")/NL*100, 1)
    rc2 = round(sum(as.numeric(tmp_text)>1)/NL*100, 1)
    bv = sum(text_left!="" & text_right!="")
    
  # Create a data frame for plotting

    cleavage_data <- data.frame(
      position = 1:NL,
      nucleotide = base::strsplit(seq1, "")[[1]],
      value = "Sequence",
      text_left = text_left,
      text_right = text_right,
      site = as.factor(cleavage_sites),
      site1 = as.factor(cleavage_sites1),
      site2 = as.factor(cleavage_sites2)
    )

  # Plot the oligonucleotide sequence and cleavage sites with nucleotide letters

    if (return.plot){
    gpl = ggplot(cleavage_data, aes(x = position, y = 0, group = value)) +
    geom_segment(data = subset(cleavage_data, site == "Cleavage Site"), aes(x = position + 0.5, y = -0.017, xend = position + 0.5, yend = 0.017), color = "red", linewidth = 1.2) +
    geom_segment(data = subset(cleavage_data, site1 == "Cleavage Site left"), aes(x = position - 0.2, y = 0.017, xend = position + 0.5, yend = 0.017),  color = "red", linewidth = 1.2) +
    geom_segment(data = subset(cleavage_data, site2 == "Cleavage Site right"), aes(x = position + 0.5, y = -0.017, xend = position + 1.2, yend = -0.017),  color = "red", linewidth = 1.2) +
    geom_text(aes(label = nucleotide), vjust = 0, size  = 8, fontface='bold') +
    geom_text(aes(x = position-0.05, y = 0.022,label = text_left), size  = 5) + 
    geom_text(aes(x = position+1.05, y = -0.020,label = text_right), size  = 5) + 
    labs(title = "McLuckey product ions", x = "", y = "", size = 10) + #ylim(-0.025, 0.03) + 
    theme_minimal() + theme(panel.grid=element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position="none", plot.title=element_text(size=30)) +   
    scale_x_continuous(expand = c(0.03, 0.03)) +  scale_y_continuous(expand = c(0.01, 0.01))}
  }
  
  coverage =  paste0("Sequence coverage is ", rc1, "% when one diagonostic fragment detected for each position and ", rc2, "% when two diagonostic fragments have to be detected. ", bv, " out of ", NL, " positions received bidirectional verification.")
  
  if (int.frag & NL>2 & nrow(features2)>0){
    
    nucleotide = base::strsplit(seq1, "")[[1]]
    kl = 1:nrow(features2)
    sites = lapply(internals, function(x) strsplit(x, "-")[[1]])
    x1 = NL- sapply(sites, function(x) as.numeric(gsub("[^0-9.-]", "", x[1])))
    x2 = sapply(sites, function(x) as.numeric(gsub("[^0-9.-]", "", x[2])))
    if (return.plot){
    gpk = ggplot() + 
      geom_segment(aes(x = x1, y = 0.023 + 0.002*kl, xend = x2, yend = 0.023 + 0.002*kl), color = "blue", linewidth = 1) +
      geom_text(aes(x = x2 + 1, y = 0.023 + 0.002*kl,label = internals), size  = 6, color = "blue") + 
      geom_text(aes(x= 1:NL, y = 0.02, label = nucleotide), vjust = 0, size  = 8, fontface='bold') +
      labs(title = "Internal fragments", x = "", y = "", size = 10) + 
      theme_minimal() + 
      theme(panel.grid=element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position="none",plot.title=element_text(size=30)) + 
      scale_x_continuous(expand = c(0.03, 0.03)) +  
      scale_y_continuous(expand = c(0.01, 0.01))}
    
  }
  
  
  return(list(gpl = gpl, gpk = gpk, cleavage_data = cleavage_data, coverage = coverage, features1 = features1, features2 = features2))
}
