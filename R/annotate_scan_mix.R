#' Targeted followed by non-targeted screening from a deconvoluted oligonucleotide spectra
#'
#' The function first searches a complex deconvoluted oligonucleotide spectra against a user provided oligonucleotide impurity/metabolite database, annotating and scoring isotopic pattern matches. It then searches DNA/RNA-like isotope patterns from the rest of deconvoluted oligonucleotide spectra. It provides the monoisotopic molecular weight, average, intensity and envelope likeness of all features detected. 
#' 
#' @author Youzhong Liu, \email{liu-youzhong@hotmail.com}
#' 
#' @param scan_processed_aggregated Data frame representing deconvoluted NMS with the true molecular weight scale. Output of the function process_scan.
#' @param MSMS Boolean. TRUE if the spectrum is MS/MS.
#' @param ntheo Integer. Estimated isotope envelop size in number of isotope peaks.
#' @param formula_flp Character. Neutral elemental formula of the main compound (e.g. full length product). Used when a transformation list from the main compound is defined. 
#' @param cpd_flp Character. Name of the main compound. Used to label compounds in the output table. Use the separator & if you expect multiple main compound in your sample.
#' @param transformation_list Data frame. Transformation list defining the mass and elemental formula difference from the main compound. Should contain following columns: ID, CPD (IDs and names of transformation products), Plus_Formula, Minus_Formula (Elemental formula difference), Delta.AVG.MW, Delta.MONO.MW (Mass difference for average and mono molecular weight)
#' @param mdb Data frame. You can directly provide all expected compounds to be annotated from your mixture without defining FLP or compound names. 
#' @param bblock Character. Either "DNA" or "RNA". Should reflect the main nucleic acid composition of the strand. Used for monoisotopic peak prediction by Pointless algorithm. 
#' @param min_overlap Double between 0 and 1. The minimum matching score between experimental and theoretical isotope envelops (known compounds from database/transformation list or unknown predicted by Pointless). 
#' @param max_msigma  Double between 1 and 50. The maximum-allowed deviation between the shapes of experimental and theoretical isotope pattern. Should set higher for noisy or MS/MS data.
#' @param max_mmw_ppm Double between 1 and 50. The maximum allowed ppm error between masses in the NMS and theoretical molecular weight of oligonucleotide features. Depend on experimental mass deviation and deconvolution bias. 
#' @param baseline Numeric. Estimated baseline level (noise) of input spectrum. Depending on instrument and acquisition method. Baseline of MS/MS spectrum is 100 for most instruments. 
#' 
#' @importFrom plyr rbind.fill
#' @importFrom BRAIN calculateAverageMass
#' @export
#'  
#' @examples
#'
#' \dontrun{ 
#' 
#' ## Example of MS1 data:
#' 
#' data("Strand_A")
#' 
#' scan.deconvoluted = process_scan(scan.A,
#' polarity = "Negative", baseline = 1000, mz_error = 0.01, 
#' min_charge = 3, max_charge = 12,
#' min_mz = 500, max_mz = 1200, min_mw = 4000, max_mw = 10000,
#' mw_gap = 1.1, mw_window = 10)
#' SCAN_NMS = scan.deconvoluted$scan_processed_aggregated
#' 
#' data("TRANS") # Transformation list for potential oligonucleotide degradants
#' 
#' scan.deconvoluted.annotated = annotate_scan_mix(SCAN_NMS, MSMS = F, ntheo = 10,
#' formula_flp = "C189H238O119N66P18S4F8",cpd_flp = "Demo A", transformation_list = transformation_list, mdb = NULL, 
#' bblock = "RNA", min_overlap = 0.6, max_msigma = 5, max_mmw_ppm = 10, baseline = 1000)
#' 
#' view(scan.deconvoluted.annotated$feature)
#' 
#' ## Example of MS2 data:
#' 
#' data("Strand_B_MS2") 
#' 
#' scan.deconvoluted = process_scan(scan2_B, MSMS = T,
#' polarity = "Negative", baseline = 100, mz_error = 0.02,
#' min_charge = 1, max_charge = 10,
#' min_mz = 500, max_mz = 1200, min_mw = 0, max_mw = 7000,
#' mw_gap = 1.1, mw_window = 6)
#' SCAN_NMS = scan.deconvoluted$scan_processed_aggregated
#' 
#' seq = "OH-Am*-Af*-Cm*-Af-Um-Uf-Gm-Af-Gm-Cf-Gm-Af-Um-Gf-Um-Cf-Cm-Am*-Cm-OH"
#' mDB = predict_esi_frag(seq) # Prediction of product ions based on sequence
#'
#'scan.deconvoluted.annotated = annotate_scan_mix(SCAN_NMS, MSMS = T, ntheo = 6,
#'          formula_flp = "",  cpd_flp = "", transformation_list = NULL, mdb = mDB, 
#'          bblock= "RNA", min_overlap = 0.4, max_msigma = 20, max_mmw_ppm = 10, baseline = 50)
#' 
#'view(scan.deconvoluted.annotated$feature)
#'
#'}
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
          feature_to_add$THEO_AMW = round(calculateAverageMass(ListFormula1(scan_to_add$FORMULA)),4)
          feature_to_add$EXP_MMW =  round(scan_to_add$MW, 4)
          feature_to_add$EXP_AMW =  -1
          feature_to_add$EXP_MMW_PPM =  round(EXP_MMW_PPM_DB[idx], 2)
          feature_to_add$EXP_AMW_DEV =  -1
          feature_to_add$SCORE = max_msigma
          feature_to_add$OC =  round(1/ntheo, 1)
          feature_to_add$RESPONSE = round(scan_to_add$Response,0)
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


#######################
##Additional functions#
#######################

ListFormula1 <- function(elemental.formula){
  chr <- gregexpr("[[:upper:]][[:lower:]]{0,1}", elemental.formula)
  for (i in 1:length(chr[[1]])) {
    y <- attr(chr[[1]], which = "match.length")[i]
    z <- substr(elemental.formula, chr[[1]][i], chr[[1]][i] +
                  y - 1)
    if (!(z == "C" | z == "H" | z == "N" |
          z == "O" | z == "S" | z == "P" |
          z == "Br" | z == "Cl" | z == "F" |
          z == "I" | z == "Si" | z == "Sn" | 
          z == "B" | z == "Na" | z== "K" | z== "Fe"))
      stop(paste("Elemental formula", elemental.formula,
                 "contains element not of C,H,N,O,S,P,Br,Cl,F,I,Si,Sn,B,Fe, Na,K."))
  }
  GetAtoms <- function(elemental.formula, element) {
    reg.exp <- paste(element, "[[:digit:]]*(?![[:lower:]])",
                     sep = "")
    #reg.exp <- paste(element, "[[:digit:]]+\\.*[[:digit:]]*(?![[:lower:]])",
    #                 sep = "")
    x <- gregexpr(reg.exp, elemental.formula, perl = TRUE)
    if (x[[1]][1] != -1) {
      n <- vector(mode = "numeric", length = length(x[[1]]))
      for (i in 1:length(x[[1]])) {
        y <- attr(x[[1]], which = "match.length")[i]
        z <- substr(elemental.formula, x[[1]][i], x[[1]][i] +
                      y - 1)
        number <- as.numeric(strsplit(z, split = element)[[1]][2])
        if (is.na(number)) {
          n[i] <- 1
        }
        else {
          n[i] <- number
        }
        atoms <- sum(n)
      }
    }
    else {
      atoms <- 0
    }
    return(atoms)
  }
  elements <- c("C", "H", "N", "O",
                "S", "P", "Br", "Cl", "F",
                "I", "Si", "Sn", "B", "Na", "K", "Fe")
  result <- as.list(sapply(elements, function(x) {
    GetAtoms(elemental.formula, x)
  }))
  return(result)
}
