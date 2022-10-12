#' Non-targeted screening from a deconvoluted oligonucleotide spectra
#'
#' The function searches DNA/RNA-like isotope patterns from a deconvoluted oligonucleotide spectra. It provides the monoisotopic molecular weight, average, intensity and envelope likeness of all features detected. 
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @export
#' 

annotate_scan_untargeted<-function(scan_processed_aggregated, bblock, ntheo = 12, 
min_relative = 0.1, min_overlap = 0.6, max_msigma = 1){
  
  # Annotate envelop
  
  features_annotated = c()
  scan_annotated = c()
  
  if (!is.null(scan_processed_aggregated)){
    EEE = max(scan_processed_aggregated$Envelop)
    for (i in 1:EEE){
      inds = which(scan_processed_aggregated$Envelop==i)
      if (length(inds)>=2){ # At least 2 peaks in the envelop
        envelop = scan_processed_aggregated[inds,]
        results = annotate_envelop_bb(envelop, bblock, ntheo, min_relative, min_overlap, max_msigma)
        features_annotated = rbind.data.frame(features_annotated, results$features.annotated)
        scan_annotated = rbind.data.frame(scan_annotated,results$envelop.annotated)
      }
    }
    
    if (!is.null(features_annotated)){
      features_annotated$AMW_EXP = round(as.numeric(features_annotated$AMW_EXP), 3)
      features_annotated$AMW_THEO = round(as.numeric(features_annotated$AMW_THEO), 3)
      features_annotated$SCORE = round(as.numeric(features_annotated$SCORE), 2)
      
    }}
  
  return(list(scan = scan_annotated, feature = features_annotated))
}


###################################################################
####Annotate envelops by pattern matching from a building block ###
###################################################################

annotate_envelop_bb<-function(envelop, bblock, ntheo = 12, 
          min_relative = 0.01, min_overlap = 0.6, max_msigma = 1){
  
  # Filter envelop:
  
  tmp_scan = cbind.data.frame(envelop$MW, envelop$Response)
  tmp_scan1 = tmp_scan
  colnames(tmp_scan1) = c("mz", "intensity")
  tmp_scan1[,2] = tmp_scan1[,2]/max(tmp_scan1[,2])*100
  valid = which(tmp_scan1[,2]>min_relative)
  tmp_scan1 = tmp_scan1[valid,]
  tmp_scan = tmp_scan[valid,]
  NP = nrow(tmp_scan)
  tic0 = sum(tmp_scan[,2]) # All scan intensity
  
  envelop = envelop[valid,]
  envelop_annotated = envelop
  
  mSigma = NULL
  
  # Decomposition for all starting MWs
  
  if (NP<=1){
    return(NULL)
  } else if (NP==2) {
    dm = tmp_scan[2,1] - tmp_scan[1,1]
    dev = tmp_scan1[2,2] - tmp_scan1[1,2]
    if (abs(dm - 1)<0.1 & dev<0){
      tmp_feature_mw = tmp_scan1[1,1]
      tmp_feature_response = sum(tmp_scan[,2])
      tmp_feature_all_mw = paste0(round(tmp_scan1[,1],2), collapse = ":")
      tmp_avg_mw = sum(tmp_scan[,1] * tmp_scan[,2]/sum(tmp_scan[,2]))
      
      envelop_annotated$SCORE = 1
      envelop_annotated$OVERLAP = 0.5
      envelop_annotated$MMW.Starter = c(1,0)
      
      sd1 = cbind.data.frame(FEATURE = paste0("E",envelop$Envelop[1]), MMW = tmp_feature_mw,
            AMW_EXP = tmp_avg_mw, AMW_THEO = tmp_avg_mw, SCORE = 1.0, RESPONSE = tmp_feature_response, MW_ALL = tmp_feature_all_mw)
      envelop_annotated$FEATURE = sd1$FEATURE
      
      list(envelop.annotated = envelop_annotated, features.annotated = sd1)      
    } else {return(NULL)}
    
  }  else {
    
    start1 = 1:max(1, NP-2)
    NS = length(start1)
    mw0_s1 = tmp_scan1[start1,1]
    mSigma_list = overlap_list = rep(0, NS)
    theo_list1 = list()
    
    for (s in start1){
      
      MW01 = mw0_s1[s]
      if (bblock %in% c("DNA", "RNA")){
        theo_isotope1 <- brain_from_pointless(type = bblock, MW01, ntheo)
      } else {
        theo_isotope1 <- brain_from_bblock(bblock, MW01, ntheo)
      }
      
      # Evaluate 
      
      mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope1, ntheo)
      ip1 = match(round(theo_isotope1$masses, 1), round(tmp_scan1[,1], 1))
      ip1 = ip1[!is.na(ip1)]
      tmp_scan2 = tmp_scan[tmp_scan[,1]>=min(tmp_scan[ip1,1]) &  tmp_scan[,1]<=max(tmp_scan[ip1,1]),]
      oc1 = sum(tmp_scan[ip1,2])/sum(tmp_scan2[,2]) # explained intensity
      
      ip2 = match(round(tmp_scan1[,1], 1), round(theo_isotope1$masses, 1))
      ip2 = ip2[!is.na(ip2)]
      oc2 = sum(theo_isotope1$isoDistr[ip2]) # explained intensity
      
      oc = 2*(oc1*oc2)/(oc1+oc2)
      
      mSigma_list[s] = mSigma$score
      overlap_list[s] = round(oc,2)
      theo_list1[[s]] = theo_isotope1
    }
    
    valid = which(overlap_list>=min_overlap & mSigma_list<=max_msigma) 
    # percent intensity explained and msigma
      
    envelop_annotated$SCORE = c(mSigma_list, 100, 100)
    envelop_annotated$OVERLAP = c(overlap_list,0, 0)
    envelop_annotated$MMW.Starter = 0
    envelop_annotated$FEATURE = "N/A"
    
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
      
      for (t in 1:(length(inds)-1)){
        all_ind = inds[t]:(inds[t+1]-1)
        tmp_theo = theo_list1[[inds[t]]]
        tmp_feature_mw = tmp_scan[inds[t],1]
        tmp_feature_response  = sum(tmp_scan[all_ind, 2])
        tmp_feature_aw = sum(tmp_scan[all_ind, 1] * tmp_scan[all_ind,2]/sum(tmp_scan[all_ind,2]))
        tmp_feature_all_mw = paste0(round(tmp_scan[all_ind, 1],2), collapse = ":")
        
        tmp_sd = cbind.data.frame(FEATURE = paste0("E",envelop$Envelop[1], "T", t), MMW = tmp_feature_mw,
                  AMW_EXP = tmp_feature_aw, AMW_THEO = tmp_theo$avgMass, 
                  SCORE = mSigma_list[inds[t]], RESPONSE = tmp_feature_response, MW_ALL = tmp_feature_all_mw)
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

calcul_imp_mSigma <- function(sp_deconvoluted, theo_isotope, ntheo, unknown = F){
  
  # Evaluate isotope patterns
  # sp_deconvoluted: m/z, intensity two column matrix
  
  # compute theoretical isotope distribution for all compounds
  
  theo_deconvoluted = cbind(theo_isotope$masses,theo_isotope$isoDistr)
  theo_deconvoluted = as.data.frame(theo_deconvoluted)
  theo_deconvoluted[,2] = theo_deconvoluted[,2]/max(theo_deconvoluted[,2])
  NT = nrow(theo_deconvoluted)
  
  sp_deconvoluted = as.data.frame(sp_deconvoluted)
  NS = nrow(sp_deconvoluted)
  
  colnames(sp_deconvoluted) = colnames(theo_deconvoluted) = c("MW","I")
  
  # Generate comparison vectors:
  
  em_list = rep(1, NT) # Dalton error from matched isotope
  res_list = rep(0, NT) # Response list from matched isotope
  exp_list = rep(0, NT) # Experimental mass
  
  for (j in 1:NT){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted$MW[j])
    valid = which(errors<0.1)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list[j] = errors[valid]
      res_list[j] = sp_deconvoluted$I[valid]
      exp_list[j] = sp_deconvoluted$MW[valid]
    }
  }
  
  filter = which(res_list>0)
  em_list = em_list[filter]
  res_list = res_list[filter]
  exp_list = exp_list[filter]
  tm_list = theo_deconvoluted$MW[filter]
  theo_list = theo_deconvoluted$I[filter]
  NT = nrow(theo_deconvoluted)
  
  # Compute evaluations:
  
  res_list = res_list/sum(res_list)
  theo_list = theo_list/sum(theo_list)
  
  chi_squa_score = sum(abs(res_list-theo_list)^2/theo_list)/NT
  
  # Sum up:
  
  em_score = sum(em_list/tm_list*1000000)/NT 
  
  return(list(score = chi_squa_score + em_score,
              exp_sp = cbind(Mass = exp_list, Res = res_list)))
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
  
  param_file <- paste0("params", type, ".rds")
  if (!file.exists(param_file)) stop("The file with model parameters could not be find.")
  model <- readRDS(param_file)
  
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
    best_mw = min_mw
    
    if ((mwlist[k] - best_mw > mw_window & mwlist[k] - max_mw > 0.8) || (mwlist[k] - max_mw >=1.1)){
      mw_feature[t0:(k-1)] = f
      mw_avg[t0:(k-1)] = best_mw
      f = f + 1
      t0 = k
    }
  }
  
  ttt=  t0:N
  mw_feature[t0:N] = f
  #mw_avg[t0:N] = mwlist[ttt[which.min(scorelist[ttt])]]
  mw_avg[t0:N] = min(mwlist[ttt])
  return(list(id = mw_feature, mw = mw_avg))
}