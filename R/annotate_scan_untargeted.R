#' Non-targeted screening from a deconvoluted oligonucleotide spectra
#'
#' The function searches DNA/RNA-like isotope patterns from a deconvoluted oligonucleotide spectra. It provides the monoisotopic molecular weight, average, intensity and envelope likeness of all features detected. 
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom BRAIN useBRAIN
#' @export
#' 

annotate_scan_untargeted<-function(scan_processed_aggregated, bblock, ntheo = 12, 
          min_overlap = 0.6, max_msigma = 10, max_mmw_ppm = 10, baseline = 1000){
  
  # Annotate envelop
  
  features_annotated = c()
  scan_annotated = c()
  
  if (!is.null(scan_processed_aggregated)){
    EEE = max(scan_processed_aggregated$Envelop)
    for (i in 1:EEE){
      inds = which(scan_processed_aggregated$Envelop==i)
      if (length(inds)>=2){ # At least 2 peaks in the envelop
        envelop = scan_processed_aggregated[inds,]
        results = annotate_envelop_bb(envelop, bblock, ntheo, min_overlap, max_msigma, max_mmw_ppm, baseline)
        features_annotated = rbind.data.frame(features_annotated, results$features.annotated)
        scan_annotated = rbind.data.frame(scan_annotated,results$envelop.annotated)
      }
    }
  }
  
 #colname_to_remove = c("OVERLAP", "AVG_MASS", "AVG_MASS_REF", "AVG_MASS_DEV")
  colname_to_remove = c("AVG_MASS", "AVG_MASS_REF", "AVG_MASS_DEV")
  inds = which(!(colnames(scan_annotated) %in% colname_to_remove)) 
  scan_annotated = scan_annotated[,inds]
  
  return(list(scan = scan_annotated, feature = features_annotated))
}

###################################################################
####Annotate envelops by pattern matching from a building block ###
###################################################################

annotate_envelop_bb<-function(envelop, bblock, ntheo = 12, min_overlap = 0.6, max_msigma = 5, max_mmw_ppm=10, baseline =  1000){
  
  # Filter envelop:
  
  tmp_scan = cbind.data.frame(envelop$MW, envelop$Response)
  tmp_scan1 = tmp_scan
  colnames(tmp_scan1) = c("mz", "intensity")
  tmp_scan1[,2] = tmp_scan1[,2]/max(tmp_scan1[,2])*100
  NP = nrow(tmp_scan)
  tic0 = sum(tmp_scan[,2]) # All scan intensity
  
  envelop_annotated = envelop

  mSigma = NULL
  
  # Decomposition for all starting MWs
  
  if (NP<=1){
    return(NULL)
  } else if (NP==2) {

    # Two peak situation - check for smaller fragments:
    
    dm = tmp_scan[2,1] - tmp_scan[1,1]
    dev = tmp_scan1[2,2] - tmp_scan1[1,2]
    
    if (abs(dm - 1)<0.1 & dev<0){
      
      tmp_feature_mw = tmp_scan1[1,1]
      tmp_feature_response = sum(tmp_scan[,2])
      #tmp_feature_all_mw = paste0(round(tmp_scan1[,1],2), collapse = ":")
      tmp_avg_mw = sum(tmp_scan[,1] * tmp_scan[,2]/sum(tmp_scan[,2]))
    
      envelop_annotated$SCORE = max_msigma
      envelop_annotated$OC_SCORE = 1
      envelop_annotated$MMW.Starter = c(1,0)
      envelop_annotated$AVG_MASS = c(round(tmp_avg_mw,4), 0)
      envelop_annotated$AVG_MASS_REF = 0
      envelop_annotated$AVG_MASS_DEV = -1
      
      sd1 = cbind.data.frame(FEATURE = paste0("E",envelop$Envelop[1]), FORMULA = "Unknown",
                          THEO_MMW = -1, THEO_AMW = -1, 
                          EXP_MMW = tmp_feature_mw, EXP_AMW =  round(tmp_avg_mw, 3),
                          EXP_MMW_PPM = -1, EXP_AMW_DEV = -1,
                          SCORE = max_msigma, OC = 1, RESPONSE = tmp_feature_response, 
                          #MW_ALL = tmp_feature_all_mw, 
                          Envelop = envelop$Envelop[1])
      envelop_annotated$FEATURE = sd1$FEATURE
      return(list(envelop.annotated = envelop_annotated, features.annotated = sd1))
    } else {return(NULL)}
    
  }  else {
    
    start1 = 1:max(1, NP-2)
    NS = length(start1)
    mw0_s1 = tmp_scan1[start1,1]
    mSigma_list = rep(100, NP)
    avg_mass_list = avg_mass_ref_list = rep(0, NP)
    overlap_list = avg_mass_dev_list = rep(-1, NP)
    theo_list1 = list()
    
    # Test different starting points for mono-isotopic!
    
    for (s in start1){
      
      MW01 = mw0_s1[s]
      if (bblock %in% c("DNA", "RNA")){
        theo_isotope1 <- brain_from_pointless1(type = bblock, MW01, ntheo)
      } else {
        theo_isotope1 <- brain_from_bblock(bblock, MW01, ntheo)
      }
      
      # Evaluate 
      
      mSigma = calcul_imp_mSigma_unknown(tmp_scan, theo_isotope1, ntheo, baseline, max_mmw_ppm)

      mSigma_list[s] = round(mSigma$score,2)
      overlap_list[s] = round(mSigma$oc_score,2)
      avg_mass_list[s] = round(mSigma$avg_mass, 4)
      avg_mass_ref_list[s] = round(mSigma$avg_mass_ref, 4)
      avg_mass_dev_list[s] = round(mSigma$avg_mass_dev,2)
      theo_list1[[s]] = theo_isotope1
    }
    
    valid = which(overlap_list>=min_overlap & mSigma_list<=max_msigma) 
    # percent intensity explained and msigma
      
    envelop_annotated$SCORE = mSigma_list
    envelop_annotated$OC_SCORE = overlap_list
    envelop_annotated$MMW.Starter = 0
    envelop_annotated$AVG_MASS = avg_mass_list
    envelop_annotated$AVG_MASS_REF = avg_mass_ref_list
    envelop_annotated$AVG_MASS_DEV = avg_mass_dev_list
  
    # Confirm monoisotopic MW: 
    
    start_mws = c()
    
    if (length(valid)==0){return(NULL)}
    if (length(valid)==1) {
      start_mws = tmp_scan[valid,1]
      envelop_annotated$MMW.Starter[valid] = 1
    } 
    if (length(valid)>1) {
      # Multiple!
      tmp_scan2 = tmp_scan[valid,,drop=FALSE]
      mSigma_list2 = mSigma_list[valid]
      envelop_annotated$MMW.Starter[valid] = 1
      start_mws = unique(cut_mmw_list3(tmp_scan2[,1], mSigma_list2, ntheo)$mw)
    }
    
    if (length(start_mws)>0){  
      
      inds = c(as.numeric(match(start_mws, tmp_scan[,1])), NP+1)
      sd2 = c()
      envelop_annotated$FEATURE = "N/A"
      
      for (t in 1:(length(inds)-1)){
        all_ind = inds[t]:(inds[t+1]-1)
        tmp_theo = theo_list1[[inds[t]]]
        tmp_feature_mw = tmp_scan[inds[t],1]
        tmp_feature_response  = sum(tmp_scan[all_ind, 2])
        tmp_feature_aw = envelop_annotated$AVG_MASS[inds[t]]
        tmp_feature_aw_ref = envelop_annotated$AVG_MASS_REF[inds[t]]
        tmp_feature_aw_dev = envelop_annotated$AVG_MASS_DEV[inds[t]]
        tmp_msigma = envelop_annotated$SCORE[inds[t]]
        tmp_oc = envelop_annotated$OC_SCORE[inds[t]]
        #tmp_feature_all_mw = paste0(round(tmp_scan[all_ind, 1],2), collapse = ":")
        
        tmp_sd = cbind.data.frame(FEATURE = paste0("E",envelop_annotated$Envelop[1], "T", t), FORMULA = "Unknown",
                          THEO_MMW = -1, THEO_AMW = tmp_feature_aw_ref,        
                          EXP_MMW = tmp_feature_mw, EXP_AMW =  round(tmp_feature_aw,3), 
                          EXP_MMW_PPM = -1, EXP_AMW_DEV = tmp_feature_aw_dev,
                          SCORE = tmp_msigma, OC = round(tmp_oc, 1), RESPONSE = tmp_feature_response, 
                          #MW_ALL = tmp_feature_all_mw,  
                          Envelop = envelop$Envelop[1])
        sd2 = rbind.data.frame(sd2, tmp_sd)
        envelop_annotated$FEATURE[all_ind] =  tmp_sd$FEATURE[1]
      }
      
      return(list(envelop.annotated = envelop_annotated, features.annotated = sd2))
    }
  }
}

###################################
### Isotope evaluation functions###
###################################

calcul_imp_mSigma_unknown <- function(sp_deconvoluted, theo_isotope, ntheo, baseline, max_mmw_ppm){
  
  # Evaluate isotope patterns
  # sp_deconvoluted: m/z, intensity two column matrix
  # Compare to targeted: overlapping score added
  
  # compute theoretical isotope distribution for all compounds
  
  theo_deconvoluted = cbind.data.frame(theo_isotope$masses,theo_isotope$isoDistr)
  theo_deconvoluted[,2] = theo_deconvoluted[,2]/sum(theo_deconvoluted[,2])
  
  tmw = theo_deconvoluted[1,1] # Theoretical mono
  taw = sum(theo_deconvoluted[,1] * theo_deconvoluted[,2])
  
  sp_deconvoluted = as.data.frame(sp_deconvoluted)
  NS = nrow(sp_deconvoluted)
  colnames(sp_deconvoluted) = colnames(theo_deconvoluted) = c("MW","I")
  NT = nrow(theo_deconvoluted)

  # Generate comparison vectors:
  
  em_list = rep(1, NT) # Dalton error from matched isotope
  res_list = rep(baseline/2, NT) # Response list from matched isotope
  exp_list = rep(0, NT) # Experimental mass
  abs_dev = max_mmw_ppm/1000000*median(sp_deconvoluted[,1])
  
  for (j in 1:NT){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted$MW[j])
    valid = which(errors<abs_dev)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list[j] = errors[valid]
      res_list[j] = sp_deconvoluted$I[valid]
      exp_list[j] = sp_deconvoluted$MW[valid]
    }
  }
  
  # Calculate chi_square, overlapping and shape similarity score:
  
  tm_list = theo_deconvoluted$MW
  theo_list = theo_deconvoluted$I
  idx = which(res_list>baseline)
  
  #oc1 = sum(res_list[idx])/sum(sp_deconvoluted[,2]) # Unnormalized
  #oc2 = sum(theo_deconvoluted$I[idx])/sum(theo_deconvoluted$I)
  #oc_score = 2*(oc1*oc2)/(oc1+oc2)
  
  oc1 = length(idx)/nrow(sp_deconvoluted)
  oc2 =  length(idx)/NT
  oc_score = 2*(oc1*oc2)/(oc1+oc2)
    
  res_list = res_list/sum(res_list)
  theo_list = theo_list/sum(theo_list)
  
  chi_squa_score = sum(abs(res_list-theo_list)/theo_list)/NT
  if (is.null(chi_squa_score)){chi_squa_score = 100}
  
  # Compute evaluations:
  
  mono_mass = exp_list[1]
  avg_mass = sum(exp_list[idx] * res_list[idx]/sum(res_list[idx])) # Average molecular weight measured
  taw_bis = sum(tm_list * theo_list/sum(theo_list)) # Theoretical average in the filtered range
  avg_mass_dev = abs(avg_mass - taw_bis) # Average mass error in absolute!!!
  
  # Sum up:
  
  if (avg_mass>0){
    return(list(score = chi_squa_score, oc_score = oc_score, 
                mono_mass_ref = -1,  mono_mass =  mono_mass,
                avg_mass_ref = taw_bis, avg_mass = avg_mass, 
                mono_ppm_dev = -1, avg_mass_dev = avg_mass_dev,
                exp_sp = cbind(Mass = exp_list, Res = res_list)))
}}

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

cut_mmw_list3<-function(mwlist, scorelist, mw_window){
  
  # More used to cluster/separate envelops!  
  
  N=length(mwlist)
  
  f=1
  mw_feature = rep(0, N)
  mw_avg = rep(0, N)
  t0 = 1 # Start index of a cluster
  
  for (k in 2:N){
    
    ttt=  t0:(k-1)
    min_mw = min(mwlist[ttt])
    max_mw = max(mwlist[ttt])
    #best_mw = mwlist[ttt[which.min(scorelist[ttt])]]
      
    if ((mwlist[k] - min_mw > 2*mw_window) || (mwlist[k] - max_mw >1.1)){ 
      mw_feature[t0:(k-1)] = f
      mw_avg[t0:(k-1)] = min_mw
      f = f + 1
      t0 = k
    }
  }
  
  ttt=  t0:N
  mw_feature[t0:N] = f
  mw_avg[t0:N] = min_mw
  return(list(id = mw_feature, mw = mw_avg))
}
