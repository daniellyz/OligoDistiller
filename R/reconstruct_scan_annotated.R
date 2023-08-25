#' Reconstruct a theoretical mass spectra based on oligonucleotide features detected
#'
#' The function provides a way to check the detected features against an original spectrum. 
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom BRAIN useBRAIN
#' @importFrom stringr str_detect str_replace_all
#' @export

reconstruct_scan_annotated<-function(scan0, scans.deconvoluted.annotated, polarity = "Negative",
      mode = c("targeted","untargeted", "mixed"), mz_error = 0.01, input_charges = 5:12, bblock = "C21 H26.4 O13.2 N7.3 P2 S0.4 F0.9", ntheo = 12){
  
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
  residual_scan = c()
  
  Features1 = features1$FEATURE
  NF =  length(Features1)
  
  for (f in 1:NF){
    
    # Theoretical distribution:
    
    if (mode == "targeted"){
      tmp_feature_label = str_replace_all(Features1[f], "[^[:alnum:]]", "")
      tmp_scan1_cpd = str_replace_all(scan1$CPD, "[^[:alnum:]]", "")
      
      valid = which(str_detect(tmp_scan1_cpd,  tmp_feature_label))
      formula1 = features1$FORMULA[f]
      theo_list = BRAIN::useBRAIN(ListFormula1(formula1), nrPeaks =  ntheo)
    } 
    if (mode == "untargeted"){
      valid = which(scan1$FEATURE == Features1[f])
      mmw =  features1$EXP_MMW[f]
      if (bblock %in% c("DNA", "RNA")){
        theo_list <- brain_from_pointless1(type = bblock, mmw, ntheo)
      } else {
        theo_list <- brain_from_bblock(bblock, mmw, ntheo)
      }
    }
    
    if (mode == "mixed"){
        tmp_feature_label = str_replace_all(Features1[f], "[^[:alnum:]]", "")
        tmp_scan1_cpd = str_replace_all(scan1$CPD, "[^[:alnum:]]", "")
        valid = which(str_detect(tmp_scan1_cpd,  tmp_feature_label))
        if (length(valid)>0){ # If is a known impurity
          formula1 = features1$FORMULA[f] 
          theo_list = BRAIN::useBRAIN(ListFormula1(formula1), nrPeaks =  ntheo)
        } else {
          valid = which(scan1$FEATURE == Features1[f])
          mmw =  features1$EXP_MMW[f]
          if (bblock %in% c("DNA", "RNA")){
            theo_list <- brain_from_pointless1(type = bblock, mmw, ntheo)
          } else {
            theo_list <- brain_from_bblock(bblock, mmw, ntheo)
          }
        }
    }
  
    theo_relative_int =  theo_list$isoDistr/sum(theo_list$isoDistr) # Normalize to sum

  # Charge summary:
    
    #ftr = scan1[valid,,drop=FALSE]
    #charges = strsplit(paste0(ftr$z, collapse = ":"),":")[[1]]
    #charges = sort(as.numeric(unique(charges)))
    charges = input_charges

  # Calculate theoretical m/z at different charge states:
      
    if (polarity == "Negative"){
      theo_mz = (theo_list$monoisotopicMass - 1.00726*charges)/charges
      theo_max = (max(theo_list$masses) - 1.00726*charges)/charges
    }
    if (polarity == "Positive"){
      theo_mz = (theo_list$monoisotopicMass + 1.00726*charges)/charges
      theo_max = (max(theo_list$masses) + 1.00726*charges)/charges
    }
        
  # Go through charge states:
      
    for (m in 1:length(theo_mz)){
      
      # Check monoisotopic peak in raw scan
      
      errors = abs(scan0[,1] - theo_mz[m])
      
      if (min(errors)<=mz_error){ 
        start_ind = which.min(errors)
        end_ind = which.min(abs(scan0[,1] - theo_max[m]))
        theo_all_int = theo_relative_int*sum(scan0[start_ind:end_ind,2]) # intensitive predicted
        raw_scan$labels[start_ind]= Features1[f]
        raw_scan$z[start_ind]= charges[m]
      
        if (polarity == "Negative"){theo_all_mz = (theo_list$masses - 1.00726*charges[m])/charges[m]}
        if (polarity == "Positive"){theo_all_mz = (theo_list$masses + 1.00726*charges[m])/charges[m]}
        
        reconstructed_mz = c(reconstructed_mz, theo_all_mz)
        reconstructed_int = c(reconstructed_int, theo_all_int)
        
        tmp_labels = tmp_z = rep("", length(theo_all_mz))
        
        tmp_labels[1] = Features1[f]
        tmp_z[1] = charges[m]
        reconstructed_label = c(reconstructed_label, tmp_labels)
        reconstructed_z = c(reconstructed_z, tmp_z)
      }
    }
  }
    
  reconstructed_scan = cbind.data.frame(mz = reconstructed_mz, int = reconstructed_int, 
                        labels = reconstructed_label, z = reconstructed_z)
  reconstructed_scan = reconstructed_scan[order(reconstructed_scan$mz),]
  
  # Find residual scan:
  
  for (r in 1:nrow(raw_scan)){
    errors = abs(as.numeric(raw_scan[r,1]) - as.numeric(reconstructed_scan[,1]))    
    if (min(errors)>mz_error){ 
        residual_scan = rbind.data.frame(residual_scan, raw_scan[r,1:2])
    }
  }

  colnames(reconstructed_scan) = colnames(raw_scan)
  colnames(residual_scan) = colnames(raw_scan)[1:2]

  return(list(original_scan = raw_scan, reconstructed_scan = reconstructed_scan,
              residual_scan = residual_scan))
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

brain_from_pointless1 <- function(type, mass, ntheo = 12) {
  
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
  if (!index) {
    #("The supplied mass falls beyond the range of the predictive model; the isotope distribution cannot be predicted.")
    output = brain_from_bblock("C29 H36 N11.5 O20.3 P3", mass, ntheo = nbPeaks)
    return(output)
  } else {
    
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
  }}


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
