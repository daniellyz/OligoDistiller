#' Targeted followed by non-targeted screening from a deconvoluted oligonucleotide spectra
#'
#' The function first searches a complex deconvoluted oligonucleotide spectra against a user provided oligonucleotide impurity/metabolite database, annotating and scoring isotopic pattern matches. It then searches DNA/RNA-like isotope patterns from the rest of deconvoluted oligonucleotide spectra. It provides the monoisotopic molecular weight, average, intensity and envelope likeness of all features detected. 
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom plyr rbind.fill
#' @export
#' 

annotate_scan_mix <- function(scan_processed_aggregated, ntheo = 10, formula_flp = "C192H239O117N73P18S4F8", cpd_flp = "Demo A",
                      transformation_list = NULL, mdb = NULL, bblock = "DNA", min_overlap = 0.6, max_msigma = 5, max_mmw_ppm = 10, baseline = 1000){
  
  options(stringsAsFactors = F)
  features.targeted = features.untargeted = features_annotated = c()
  scans.targeted = scans.untargeted = scan_anntoated = c()
  
  # Targeted screening:  

  scan.deconvoluted = scan_processed_aggregated
  scan.deconvoluted.annotated = annotate_scan_targeted(scan.deconvoluted, formula_flp = formula_flp, cpd_flp = cpd_flp, transformation_list = transformation_list, mdb = mdb, 
                                                       ntheo, min_overlap, max_msigma, max_mmw_ppm, baseline)
  features.targeted = scan.deconvoluted.annotated$feature
  scans.targeted = scan.deconvoluted.annotated$scan
  
  ind.target = which(scans.targeted$Envelop %in% features.targeted$Envelop)
  if (length(ind.target)>0){scans.targeted = scans.targeted[ind.target,,drop=FALSE]}
  
  # Unannotated envelops:  
  
  ind.rest = which(!(scan.deconvoluted$Envelop %in% features.targeted$Envelop))

  # Untargeted screening:  
  
  if (length(ind.rest)>=3){
      scan.deconvoluted.rest = scan.deconvoluted[ind.rest,,drop=FALSE]
      scan.deconvoluted.rest.anotated = annotate_scan_untargeted(
        scan.deconvoluted.rest, bblock = bblock, ntheo = ntheo, min_overlap = min_overlap, max_msigma = max_msigma, max_mmw_ppm = max_mmw_ppm, baseline = baseline)
      
      features.untargeted = scan.deconvoluted.rest.anotated$feature
      scans.untargeted = scan.deconvoluted.rest.anotated$scan
  }
  
  # Summarize:
  
  scan.total = rbind.fill(scans.targeted, scans.untargeted)
  scan.total = scan.total[order(scan.total$MW),]
  feature.total = rbind.fill(features.targeted, features.untargeted)
  feature.total = feature.total[order(feature.total$EXP_AMW),]
  
  filter = which(feature.total$SCORE<=max_msigma)
  feature.total = feature.total[filter,,drop=FALSE]
  
  return(list(scan = scan.total, feature = feature.total))
  
}


