
library(BRAIN)
library(OrgMassSpecR)
options(stringsAsFactors = F)

##############################################
####Annotate envelops by searching database###
##############################################

annotate_envelop<-function(envelop, ref_trans, IFL, ntheo = 12){
  
  # Envelop annotation:
  
  envelop_annotated = c()
  
  sd1 = envelop
  tmp_cpd = tmp_formula = tmp_score = tmp_cpd1 = rep(NA, nrow(sd1))
  tmp_feature = tmp_feature_envelop = tmp_feature_formula = tmp_feature_mw = tmp_feature_score = tmp_feature_response= c()
  
  coef1 = 1
  coef2 = 0
  
  tmp_scan = cbind.data.frame(sd1$MW, sd1$Response)
  avg_envelop = mean(sd1$MW)
  
  tmp_mw = ref_trans$AVG.MW
  avg_dev = abs(avg_envelop - tmp_mw)
  valid = which(avg_dev<=10)
  mSigma = NULL
  
  # Only 1 possibility
  
  if (length(valid)==1){
    ifl = IFL[[valid]]
    theo_isotope <- BRAIN::useBRAIN(aC = ifl, nrPeaks = ntheo)
    mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope, ntheo)
    coef1 = 1
    coef2 = 0
  }
  
  # 2 possibilities:
  
  if (length(valid)>=2){
    
    valid = sort(order(avg_dev)[1:2]) # Find top 2 match
    valid = sort(valid)
    tmp_scan1 = tmp_scan
    tmp_scan1[,1] = tmp_scan1[,1] - 1.00726
    colnames(tmp_scan1) = c("mz", "intensity")
    ifl1 = IFL[[valid[1]]]
    ifl2 = IFL[[valid[2]]]
    theo_isotope1 <- BRAIN::useBRAIN(aC = ifl1, nrPeaks = ntheo)
    theo_isotope2 <- BRAIN::useBRAIN(aC = ifl2, nrPeaks = ntheo)
    
    optim <- try(deconvolution1(scan_df = tmp_scan1,
                                theor_ID_cmpd1 = theo_isotope1,
                                theor_ID_cmpd2 = theo_isotope2,
                                n_theor_peaks = 15,
                                expected_charge_range = 1,
                                matching_mass_accuracy = 10,
                                noise_threshold = 0), silent = T)
    
    if (class(optim)!="try-error"){
      tmp_optim = "correct"
      coef2 = optim$by_charge$z1[1]
      coef1 = 1 - coef2
    } else {
      tmp_optim = as.character(strsplit(as.character(optim)[[1]], "\n")[[1]][2])
      if (tmp_optim == "  there is no overlap"){
        coef1 = 1
        coef2 = 1
      } else {
        coef1 = 1
        coef2 = 0
      }}
    
    if (coef1>0 & coef2>0){
      mSigma = calcul_mix_mSigma(tmp_scan, theo_isotope1, theo_isotope2, coef1, coef2, ntheo)
    }
    if (coef1>0 & coef2==0){
      valid = valid[1]
      mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope1, ntheo)
    }
    if (coef1==0 & coef2>0){
      valid = valid[2]
      mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope2, ntheo)
    }
  }
  
  # Output results:
  
  if (!is.null(mSigma)){
    if (!is.na(mSigma$score) & nrow(mSigma$exp_sp)>0){
     if (mSigma$score>0){
      if (!("score1" %in% names(mSigma))){
        coef1 = 1
        coef2 = 0
        inds = match(mSigma$exp_sp[,1], sd1$MW)
        inds = inds[!is.na(inds)]
        
        tmp_cpd[inds] = ref_trans$CPD[valid]
        tmp_formula[inds] = ref_trans$FORMULA[valid]
        tmp_score[inds] = round(mSigma$score,2)
        tmp_cpd1[inds] = ref_trans$CPD[valid]
        
        tmp_feature = ref_trans$CPD[valid]
        tmp_feature_envelop = sd1$Envelop[1]
        tmp_feature_formula = ref_trans$FORMULA[valid]
        tmp_feature_mw = tmp_mw[valid]
        tmp_feature_score = round(mSigma$score,2)
        tmp_feature_response = round(sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp)]),0)
      } else {
        inds1 = match(mSigma$exp_sp1[,1], sd1$MW)
        inds1 = inds1[!is.na(inds1)]
        tmp_cpd[inds1] = ref_trans$CPD[valid[1]]
        tmp_formula[inds1] = ref_trans$FORMULA[valid[1]]
        tmp_score[inds1] = round(mSigma$score1,2)
        
        inds2 = match(mSigma$exp_sp2[,1], sd1$MW)
        inds2 = inds2[!is.na(inds2)]
        tmp_cpd[inds2] = paste0(tmp_cpd[inds2], ":", ref_trans$CPD[valid[2]])
        tmp_formula[inds2] = paste0(tmp_formula[inds2], ":", ref_trans$FORMULA[valid[2]])
        tmp_score[inds2] = paste0(tmp_score[inds2], ":", round(mSigma$score2,2))
        
        inds = unique(c(inds1, inds2))
        tmp_score[inds] = paste0(tmp_score[inds], "%", round(mSigma$score,2))
        tmp_cpd1[inds] =  paste0(ref_trans$CPD[valid[1]] , ":",ref_trans$CPD[valid[2]])
        
        tmp_feature = ref_trans$CPD[valid]
        tmp_feature_envelop = sd1$Envelop[1]
        tmp_feature_formula = ref_trans$FORMULA[valid]
        tmp_feature_mw = tmp_mw[valid]
        tmp_feature_score = round(c(mSigma$score, mSigma$score),2)
        if (tmp_optim == "  there is no overlap"){
          tmp_r1 = sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp1)])
          tmp_r2 = sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp2)])
        } else {
          tmp_r1 = coef1*sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp)])
          tmp_r2 = coef2*sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp)])
        }
        tmp_feature_response = round(c(tmp_r1, tmp_r2), 0)
      }}}}
  
  sd1$CPD = as.character(tmp_cpd)
  sd1$FORMULA = as.character(tmp_formula)
  sd1$SCORE = as.character(tmp_score)
  sd1$CPD = str_remove(sd1$CPD,"NA:")
  
  sd1$FORMULA = str_remove(sd1$FORMULA,"NA:")
  sd1$SCORE = str_remove(sd1$SCORE,"NA:")
  sd1$CPD1 = as.character(tmp_cpd1)
  sd1$COEF = paste0(round(coef1,2), ":", round(coef2,2))
  
  sd2 =  cbind.data.frame(FEATURE = tmp_feature, FORMULA = tmp_feature_formula,
          AMW = tmp_feature_mw, SCORE = tmp_feature_score, RESPONSE = tmp_feature_response,
          Envelop = tmp_feature_envelop)

  return(list(envelop.annotated = sd1, features.annotated = sd2))
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

calcul_mix_mSigma <- function(sp_deconvoluted, theo_isotope1, theo_isotope2, coef1, coef2, ntheo){
  
  # Evaluate isotope patterns
  # dat: m/z, intensity two column matrix
  
  # compute and weigh theoretical isotope distribution for two compounds
  
  #theo_isotope1 <- BRAIN::useBRAIN(aC = ifl1, nrPeaks = ntheo)
  #theo_isotope2 <- BRAIN::useBRAIN(aC = ifl2, nrPeaks = ntheo)
  
  theo_deconvoluted1 = cbind.data.frame(Mass = theo_isotope1$masses,Intensity = theo_isotope1$isoDistr)
  theo_deconvoluted2 = cbind.data.frame(Mass = theo_isotope2$masses,Intensity = theo_isotope2$isoDistr)
  theo_deconvoluted1[,2] = theo_deconvoluted1[,2]/max(theo_deconvoluted1[,2])*coef1
  theo_deconvoluted2[,2] = theo_deconvoluted2[,2]/max(theo_deconvoluted2[,2])*coef2
  colnames(sp_deconvoluted) = colnames(theo_deconvoluted1) = colnames(theo_deconvoluted2) = c("MW","I")
  
  # Mix two distributions, Extract theoretical distributions from mix:
  
  dat_theo = rbind.data.frame(theo_deconvoluted1,theo_deconvoluted2)
  dat_theo = dat_theo[order(dat_theo[,1]),]
  dat_theo$group = cut_mmw_list(dat_theo$MW, dat_theo$I, ntheo)$id
  theo_mass= aggregate(dat_theo$MW, list(dat_theo$group), FUN=mean) 
  theo_int= aggregate(dat_theo$I, list(dat_theo$group), FUN=sum) 
  theo_deconvoluted = cbind.data.frame(MW = theo_mass$x, I = theo_int$x)
  theo_deconvoluted1[,2] = theo_deconvoluted1[,2]/max(theo_deconvoluted[,2])
  theo_deconvoluted2[,2] = theo_deconvoluted2[,2]/max(theo_deconvoluted[,2])
  
  # Generate comparison vectors:
  
  NT1 = nrow(theo_deconvoluted1)
  NT2 = nrow(theo_deconvoluted2)
  NT = nrow(theo_deconvoluted)
  
  em_list1 = rep(1, NT1) # Dalton error from matched isotope
  em_list2 = rep(1, NT2)
  em_list = rep(1, NT)
  res_list1 = exp_list1 = rep(0, NT1) # Response list from matched isotope
  res_list2 = exp_list2 = rep(0, NT2) # Experimental mass
  res_list = exp_list = rep(0, NT)
  
  for (j in 1:NT1){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted1$MW[j])
    valid = which(errors<0.1)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list1[j] = errors[valid]
      res_list1[j] = sp_deconvoluted$I[valid]
      exp_list1[j] = sp_deconvoluted$MW[valid]}
  }
  
  for (j in 1:NT2){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted2$MW[j])
    valid = which(errors<0.1)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list2[j] = errors[valid]
      res_list2[j] = sp_deconvoluted$I[valid]
      exp_list2[j] = sp_deconvoluted$MW[valid]}
  }
  
  for (j in 1:NT){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted$MW[j])
    valid = which(errors<0.1)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list[j] = errors[valid]
      res_list[j] = sp_deconvoluted$I[valid]
      exp_list[j] = sp_deconvoluted$MW[valid]}
  }
  

  filter1 = which(res_list1>0 & theo_deconvoluted1$I>0)
  em_list1 = em_list1[filter1]
  res_list1 = res_list1[filter1]
  exp_list1 = exp_list1[filter1]
  tm_list1 = theo_deconvoluted1$MW[filter1]
  theo_list1 = theo_deconvoluted1$I[filter1]
  
  filter2 = which(res_list2>0 & theo_deconvoluted2$I>0)
  em_list2 = em_list2[filter2]
  res_list2 = res_list2[filter2]
  exp_list2 = exp_list2[filter2]
  tm_list2 = theo_deconvoluted1$MW[filter2]
  theo_list2 = theo_deconvoluted1$I[filter2]
  
  filter = which(res_list>0 & theo_deconvoluted$I>0)
  em_list = em_list[filter]
  res_list = res_list[filter]
  exp_list = exp_list[filter]
  tm_list = theo_deconvoluted$MW[filter]
  theo_list = theo_deconvoluted$I[filter]
  
  # Compute evaluations:
  
  NT1 = length(theo_list1)
  NT2 = length(theo_list2)
  NT = length(theo_list)
  
  res_list1 = res_list1/sum(res_list1)
  res_list2 = res_list2/sum(res_list2)
  res_list = res_list/sum(res_list)
  
  theo_list1 = theo_list1/sum(theo_list1)
  theo_list2 = theo_list2/sum(theo_list2)
  theo_list = theo_list/sum(theo_list)
  
  chi_squa_score1 = sum(abs(res_list1-theo_list1)/theo_list1)/NT1
  chi_squa_score2 = sum(abs(res_list2-theo_list2)/theo_list2)/NT2
  chi_squa_score = sum(abs(res_list-theo_list)/theo_list)/NT
  
  em_list_score1 = sum(em_list1/tm_list1*1000000)/NT1
  em_list_score2 = sum(em_list2/tm_list2*1000000)/NT2
  em_list_score = sum(em_list/tm_list*1000000)/NT

  # Sum up:
  
  return(list(score1 = chi_squa_score1 + em_list_score1,
              score2 = chi_squa_score2 + em_list_score2,
              score = chi_squa_score + em_list_score,
              exp_sp1 = cbind(Mass = exp_list1, Res = res_list1),
              exp_sp2 = cbind(Mass = exp_list2, Res = res_list2),
              exp_sp = cbind(Mass = exp_list, Res = res_list)))
}

#######################
##Additional functions#
#######################

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

calcul_imp_formula<-function(formula_flp, ref_trans){
  
  NT = nrow(ref_trans)
  ifl_list = list()
  
  if (formula_flp!="N/A"){
    rfl  = unlist(ListFormula1(formula_flp))
  
    for (i in 1:NT){
      plus_formula = ref_trans$Plus_Formula[i]
      minus_formula = ref_trans$Minus_Formula[i]
    
      pfl = mfl = rep(0, 16)
      if (plus_formula!="N/A"){pfl = unlist(ListFormula1(plus_formula))}
      if (minus_formula!="N/A"){mfl = unlist(ListFormula1(minus_formula))}
    
    ifl = rfl+pfl-mfl 
    ifl_list[[i]] = as.list(ifl)
  }} else {
    for (i in 1:NT){
      ifl  = unlist(ListFormula1(ref_trans$FORMULA[i]))
      ifl_list[[i]] = as.list(ifl)
    }
  }
  
  return(ifl_list)
}

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

cut_mmw_list<-function(mwlist, intlist, mw_window){
  
  # More used to cluster individual molecular weights!  
  
  N=length(mwlist)
  
  f=1
  mw_feature = rep(0, N)
  mw_avg = rep(0, N)
  t0 = 1 # Start index of a cluster
  
  for (k in 2:N){
    
    ttt=  t0:(k-1)
    min_mw = min(mwlist[ttt])
    best_mw = mwlist[ttt[which.max(intlist[ttt])]]
    max_mw = max(mwlist[ttt])
    
    if (mwlist[k] - min_mw > mw_window || mwlist[k] - max_mw > 0.2){
      mw_feature[t0:(k-1)] = f
      mw_avg[t0:(k-1)] = round(best_mw, 4)
      f = f + 1
      t0 = k
    }
  }
  ttt = t0:N
  mw_feature[ttt] = f
  mw_avg[ttt] = round(mwlist[ttt[which.max(intlist[ttt])]],4)
  return(list(id = mw_feature, mw = mw_avg))
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


expand_transformation_list<-function(formula_flp, transformation_list){

  #transformation_list = read.csv("Transformation_list_jennifer.txt", sep = "\t", stringsAsFactors = F)
  # Oxaluric Acid (-20), Unlock nucleic acid, Deamination (+1Da), 2' Fluorination (+19)
  # Spiroiminodihydantion (+32), (S)-constrained ethyl (+26), 2'-Methoxy (O-methyl)(+31)
  # 2'-O-ally(40), Acetylation of Cytosine (+42), Locked nucleic acid (44),ADMF/GDMF (57), 2' Methoxyethyl(58)

  
  amwFLP = BRAIN::calculateAverageMass(ListFormula1(formula_flp))
  mmwFLP = BRAIN::calculateMonoisotopicMass(ListFormula1(formula_flp))
  
  IFL = calcul_imp_formula(formula_flp, transformation_list)
  NIF = length(IFL)
  scan_annotated = c()
  features_annotated = c()

  for (i in 1:length(IFL)){
    ifl = unlist(IFL[[i]]) 
    ifl = ifl[which(ifl>0)]
    ifl[ifl==1] = ""
    ifl = paste0(paste0(names(ifl), ifl),collapse = "")
    transformation_list$FORMULA[i] = ifl
  }

  transformation_list$AVG.MW = sapply(IFL, BRAIN::calculateAverageMass)
  transformation_list$MONO.MW = sapply(IFL, BRAIN::calculateMonoisotopicMass)
  transformation_list$Delta.AVG.MW = transformation_list$AVG.MW - amwFLP
  transformation_list$Delta.MONO.MW = transformation_list$MONO.MW - mmwFLP
  
  return(transformation_list)
}
