#' Targeted followed by non-targeted screening from a deconvoluted oligonucleotide spectra
#'
#' The function first searches a complex deconvoluted oligonucleotide spectra against a user provided oligonucleotide impurity/metabolite database, annotating and scoring isotopic pattern matches. It then searches DNA/RNA-like isotope patterns from the rest of deconvoluted oligonucleotide spectra. It provides the monoisotopic molecular weight, average, intensity and envelope likeness of all features detected. 
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom plyr rbind.fill
#' @export
#' 

annotate_scan_mix <- function(scan_processed_aggregated, MSMS = F, ntheo = 10, formula_flp = "C192H239O117N73P18S4F8", cpd_flp = "Demo A",
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
  feature.total = rbind.fill(features.targeted, features.untargeted)
  
  # Extra for MS/MS searching in database
  
  ind.rest = which(!(scan.deconvoluted$Envelop %in% feature.total$Envelop))
  
  if (length(ind.rest)>=3 & MSMS & !is.null(mdb)){
      
      scan.deconvoluted.rest = scan.deconvoluted[ind.rest,,drop=FALSE]
      
      for (i in 1:nrow(scan.deconvoluted.rest)){
  
        dev = scan.deconvoluted.rest$MW[i] - mdb$NM
        res = scan.deconvoluted.rest$Response[i]
        rmw = round(scan.deconvoluted.rest$MW[i], 0)
        
        EXP_MMW_PPM_DB = abs(dev)/mdb$NM*1000000
        
        if (min(EXP_MMW_PPM_DB) < max_mmw_ppm & res > 2.5*baseline){ # at least 2.5 times the baseline
          
          idx = which.min(EXP_MMW_PPM_DB)
          scan_to_add = scan.deconvoluted.rest[i,,drop=FALSE]
          scan_to_add$CPD = mdb$CPD[idx]
          scan_to_add$FORMULA = mdb$FORMULA[idx]
          scan_to_add$SCORE = max_msigma
          scan_to_add$OC_SCORE = round(1/ntheo, 1)
          scan_to_add$CPD1 =  scan_to_add$CPD
          scan_to_add$COEF = "1:0"
          scan_to_add$MMW.Starter = 1
          scan_to_add$FEATURE = paste0("E", scan_to_add$Envelop[i], "_x_", rmw)
          
          feature_to_add = feature.total[1,,drop= FALSE]
          feature_to_add$FEATURE = mdb$CPD[idx]
          feature_to_add$FORMULA = scan_to_add$FORMULA 
          feature_to_add$THEO_MMW = round(mdb$NM[idx], 4)
          feature_to_add$THEO_AMW = -1
          feature_to_add$EXP_MMW =  round(scan_to_add$MW, 4)
          feature_to_add$EXP_AMW =  -1
          feature_to_add$EXP_MMW_PPM =  round(EXP_MMW_PPM_DB[idx], 2)
          feature_to_add$EXP_AMW_DEV =  -1
          feature_to_add$SCORE = max_msigma
          feature_to_add$OC =  round(1/ntheo, 1)
          feature_to_add$RESPONSE = scan_to_add$Response
          feature_to_add$Envelop = scan_to_add$Envelop
          
          scan.total = rbind.fill(scan.total, scan_to_add)
          feature.total = rbind.fill(feature.total, feature_to_add)
        }
      }
  }
  
  # Final ordering:
  
  scan.total = scan.total[order(scan.total$MW),]
  feature.total = feature.total[order(feature.total$EXP_AMW),]
  
  filter = which(feature.total$SCORE<=max_msigma)
  feature.total = feature.total[filter,,drop=FALSE]
  
  return(list(scan = scan.total, feature = feature.total))
  
}


