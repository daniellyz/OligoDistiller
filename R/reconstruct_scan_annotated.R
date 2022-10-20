#' Reconstruct a theoretical mass spectra based on oligonucleotide features detected
#'
#' The function provides a way to check the detected features against an original spectrum. 
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom BRAIN useBRAIN
#' @export
#' 
reconstruct_scan_annotated<-function(scan0, scans.deconvoluted.annotated, polarity = "Negative", baseline = 100,
                  mode = c("targeted","untargeted"), mz_error = 0.01, bblock = "C21 H26.4 O13.2 N7.3 P2 S0.4 F0.9", ntheo = 12){
  
  raw_scan = data.frame(scan0)
  raw_scan$labels = ""
  raw_scan$z = ""
  colnames(raw_scan)[1:2] = c("m/z", "I")
  
  scan1 = scans.deconvoluted.annotated$scan
  features1 = scans.deconvoluted.annotated$feature
  reconstructed_mz = c()
  reconstructed_int = c()
  reconstructed_label = c()
  reconstructed_z = c()
  
  if (mode == "targeted"){
    envelops1 = unique(features1$Envelop)
    NE =  length(envelops1)
    
    for (e in 1:NE){
      valid = which(scan1$Envelop == envelops1[e])
      envelop = scan1[valid,,drop=FALSE]
      charges = strsplit(paste0(envelop$z, collapse = ":"),":")[[1]]
      charges = sort(as.numeric(unique(charges)))
      formulas =  unlist(sapply(envelop$FORMULA, function(x) strsplit(x, ":")[[1]]))
      formulas = unique(as.character(formulas))
      formulas = formulas[!is.na(formulas)]
      labels =  unlist(sapply(envelop$CPD, function(x) strsplit(x, ":")[[1]]))
      labels = unique(as.character(labels))
      labels = labels[!is.na(labels)]

      for (j in 1:length(formulas)){
        theo_list = BRAIN::useBRAIN(ListFormula1(formulas[j]))
        theo_relative_int =  theo_list$isoDistr/theo_list$isoDistr[1]
        
        # Calculate m/z at different charge states:
        if (polarity == "Negative"){
          theo_mz = (theo_list$monoisotopicMass - 1.00726*charges)/charges
        }
        if (polarity == "Positive"){
          theo_mz = (theo_list$monoisotopicMass + 1.00726*charges)/charges
        }
        
        # Go through charge states:
        
        for (m in 1:length(theo_mz)){
          errors = abs(scan0[,1] - theo_mz[m])
          # Use monoisotopic peak found in raw scan to calculate intensities
          if (min(errors)<=mz_error){ 
            start_ind = which.min(errors)
            theo_all_int = theo_relative_int*scan0[start_ind,2]
            raw_scan$labels[start_ind]= labels[j]
            raw_scan$z[start_ind]= charges[m]
          # Use baseline if monoisotopic m/z not found in raw data
          } else {
            theo_all_int = theo_relative_int*baseline
          }
          if (polarity == "Negative"){
            theo_all_mz = (theo_list$masses - 1.00726*charges[m])/charges[m]
          }
          if (polarity == "Positive"){
            theo_all_mz = (theo_list$masses + 1.00726*charges[m])/charges[m]
          }
          reconstructed_mz = c(reconstructed_mz, theo_all_mz)
          reconstructed_int = c(reconstructed_int, theo_all_int)
          tmp_labels = tmp_z = rep("", length(theo_all_mz))
          tmp_labels[1] = labels[j]
          tmp_z[1] = charges[m]
          reconstructed_label = c(reconstructed_label, tmp_labels)
          reconstructed_z = c(reconstructed_z, tmp_z)
        }
      }
    }
  
    reconstructed_scan = cbind.data.frame(mz = reconstructed_mz, int = reconstructed_int, 
                                        labels = reconstructed_label, z = reconstructed_z)
    reconstructed_scan = reconstructed_scan[order(reconstructed_scan$mz),]
    reconstructed_scan = reconstructed_scan[reconstructed_scan$int>=baseline,]
    colnames(reconstructed_scan) = colnames(raw_scan)
    
  }
  
  if (mode == "untargeted"){
    Features1 = features1$FEATURE
    NF =  length(Features1)
    
    for (f in 1:NF){
      valid = which(scan1$FEATURE == Features1[f])
      ftr = scan1[valid,,drop=FALSE]
      charges = strsplit(paste0(ftr$z, collapse = ":"),":")[[1]]
      charges = sort(as.numeric(unique(charges)))
      mmw =  features1$MMW[f]
      
      if (bblock %in% c("DNA", "RNA")){
        theo_list <- brain_from_pointless(type = bblock, mmw, ntheo)
      } else {
        theo_list <- brain_from_bblock(bblock, mmw, ntheo)
      }
      
      theo_relative_int =  theo_list$isoDistr/theo_list$isoDistr[1]
      
      # Calculate m/z at different charge states:
      if (polarity == "Negative"){
          theo_mz = (theo_list$monoisotopicMass - 1.00726*charges)/charges
      }
      if (polarity == "Positive"){
          theo_mz = (theo_list$monoisotopicMass + 1.00726*charges)/charges
      }
        
      # Go through charge states:
      
      for (m in 1:length(theo_mz)){
        errors = abs(scan0[,1] - theo_mz[m])
        # Use monoisotopic peak found in raw scan to calculate intensities
        if (min(errors)<=mz_error){ 
          start_ind = which.min(errors)
          theo_all_int = theo_relative_int*scan0[start_ind,2]
          raw_scan$labels[start_ind]= Features1[f]
          raw_scan$z[start_ind]= charges[m]
            # Use baseline if monoisotopic m/z not found in raw data
        } else {
          theo_all_int = theo_relative_int*baseline
        }
        if (polarity == "Negative"){
          theo_all_mz = (theo_list$masses - 1.00726*charges[m])/charges[m]
        }
        if (polarity == "Positive"){
          theo_all_mz = (theo_list$masses + 1.00726*charges[m])/charges[m]
        }
        reconstructed_mz = c(reconstructed_mz, theo_all_mz)
        reconstructed_int = c(reconstructed_int, theo_all_int)
        tmp_labels = tmp_z = rep("", length(theo_all_mz))
        tmp_labels[1] = Features1[f]
        tmp_z[1] = charges[m]
        reconstructed_label = c(reconstructed_label, tmp_labels)
        reconstructed_z = c(reconstructed_z, tmp_z)
        }
      }
    
    reconstructed_scan = cbind.data.frame(mz = reconstructed_mz, int = reconstructed_int, 
                        labels = reconstructed_label, z = reconstructed_z)
    reconstructed_scan = reconstructed_scan[order(reconstructed_scan$mz),]
    reconstructed_scan = reconstructed_scan[reconstructed_scan$int>=baseline,]
    colnames(reconstructed_scan) = colnames(raw_scan)
  }

  return(list(original_scan = raw_scan, reconstructed_scan = reconstructed_scan))
}

####################################################
### Estimate isotope shape from a building block ###
####################################################

brain_from_bblock<-function(bblock, MW0, ntheo = 12){
  
  # Estimate a approximate brain output from a molecular weight using building block
  
  element = c(12, 1.007825, 14.003, 15.9949, 31.9721, 30.9738, 78.9183, 34.9689, 18.9984, 126.9045, 27.9769,  119.9022,  11.0093,  22.9898,  38.9637, 55.9349)
  
  bblock0_v = unlist(ListFormula2(bblock))
  mw0 = sum(element * bblock0_v)
  N = MW0/mw0
  BBLOCK1 = as.list(round(N*bblock0_v))
  theor_ID_cmpd1 <- BRAIN::useBRAIN(aC = BBLOCK1, nrPeaks = ntheo)
  
  # Correct based on mw0:
  
  delta_mw = MW0 - theor_ID_cmpd1$monoisotopicMass
  theor_ID_cmpd1$monoisotopicMass = MW0
  theor_ID_cmpd1$masses = delta_mw + theor_ID_cmpd1$masses
  theor_ID_cmpd1$avgMass = delta_mw + theor_ID_cmpd1$avgMass
  
  return(theor_ID_cmpd1)
}

######################################################
### Estimate isotope shape using pointness for DNA ###
######################################################

brain_from_pointless <- function(type, mass, ntheo = 12) {
  
  nbPeaks = ntheo
  
  if (any(!is.numeric(mass), length(mass) > 1)) stop("The supplied mass vector should be a scalar")
  L <- length(mass)
  
  if (type == "DNA"){
    data(DNA)
    model = DNA
  }
  if (type == "RNA"){
    data(RNA)
    model = RNA
  }
  
  index <- model$range[1] <= mass & mass <= model$range[2]
  if (!index) stop("The supplied mass falls beyond the range of the predictive model; the isotope distribution cannot be predicted.")
  
  nbIso <- dim(model$poly_coefs)[2]
  nbPoly <- dim(model$poly_coefs)[1] - 1
  
  if (nbIso < nbPeaks | nbPeaks < 1) {
    warning(paste0("Maximum number of isotope peaks is limitted to ", nbIso, "and should be larger than 0. We proceed with the default value."))
    nbPeaks <- nbIso
  }
  # standardize mass covariate
  z <- (mass - model$mu) / model$sigma
  
  # predict ALR
  mass_mat <- matrix(rep(z, each = nbPoly + 1), ncol = nbPoly + 1, byrow = TRUE)
  mass_mat <- t(apply(mass_mat, 1, function(x) "^"(x, 0:nbPoly)))
  q_ALR_hat <- mass_mat %*% model$poly_coefs
  q_ALR_hat <- unname(q_ALR_hat)
  
  # perform anti-ALR via softmax with one added to vector
  tmp <- cbind(1, exp(q_ALR_hat))
  rs <- rowSums(tmp)
  q_hat <- t(apply(cbind(rs, tmp), 1, function(x) x[-1] / x[1]))
  
  # compute masses
  m_hat <- mass + matrix(rep(model$mass_shift, L), nrow = L, ncol = nbIso, byrow = TRUE)
  
  # avg_mass
  avg_mass <- sum(m_hat[1,] * q_hat[1,-ncol(q_hat)])/sum(q_hat[1,-ncol(q_hat)])
  
  # restrict the result to the specified number of isotope peaks
  m_hat <- m_hat[, 1:nbPeaks]
  q_hat <- q_hat[, 1:nbPeaks]
  
  return(list(isoDistr = q_hat, 
              masses = m_hat, 
              monoisotopicMass = m_hat[1],
              avgMass = avg_mass))
}

#######################
##Additional functions#
#######################

ListFormula2 <- function(elemental.formula){
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
    #reg.exp <- paste(element, "[[:digit:]]*(?![[:lower:]])",
    #                 sep = "")
    reg.exp <- paste(element, "[[:digit:]]+\\.*[[:digit:]]*(?![[:lower:]])",
                     sep = "")
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
