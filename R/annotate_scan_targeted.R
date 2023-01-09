#' Annotating a deconvoluted oligonucleotide spectra
#'
#' The function searches a complex deconvoluted oligonucleotide spectra against a user provided oligonucleotide impurity/metabolite database, annotating and scoring isotopic pattern matches. Alternatively, user can provide the FLP (full length product) formula and a list of bio or chemical transformations.
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom BRAIN calculateAverageMass calculateMonoisotopicMass useBRAIN
#' @importFrom OrgMassSpecR MolecularWeight
#' @importFrom stringr str_remove str_trim
#' 
#' @export

annotate_scan_targeted<-function(scan_processed_aggregated, formula_flp = "C192H239O117N73P18S4F8", cpd_flp = "Demo A", transformation_list = NULL, mdb = NULL, 
                ntheo = 12, min_overlap = 0.6, max_msigma = 3, max_mmw_ppm = 10, baseline = 1000){
  
  options(stringsAsFactors = F)
  features_annotated = c()
  scan_annotated = c()
  
  # Calculate impurity formula:
  
  if (formula_flp!= "" & !is.null(transformation_list)){
    formula_flp = str_trim(strsplit(formula_flp, "\\&|\\%|:")[[1]])
    cpd_flp = str_trim(strsplit(cpd_flp, "\\&|\\%|:")[[1]], side = "right")
    cpd_flp = str_trim(cpd_flp, side = "left")
    if (formula_flp[1]!="" & length(formula_flp)!= length(cpd_flp)){stop("Number of Formulas should be the same as number of strand names!")}
  
    NFF = length(formula_flp)
    if (NFF>0){
      transformation_list0 = transformation_list
      transformation_list = c()
      IFL = list()
      for (k in 1:NFF){
        expanded = expand_transformation_list(formula_flp[k], transformation_list0, mdb = NULL)  
        tmp_transformation = expanded$transformation_list
        tmp_transformation$CPD = paste0(paste0(cpd_flp[k], ": "), tmp_transformation$CPD)
        tmp_ifl = expanded$IFL
        transformation_list = rbind.data.frame(transformation_list, tmp_transformation)
        IFL = c(IFL, tmp_ifl)
    }
      transformation_list$ID = 1:nrow(transformation_list)
    }} else {
      expanded = expand_transformation_list(formula_flp, transformation_list = NULL, mdb)  
      transformation_list = expanded$transformation_list
      IFL = expanded$IFL
    }
  
  # Annotate envelop
  
  if (!is.null(scan_processed_aggregated) & nrow(transformation_list)>0){
    EEE = max(as.numeric(scan_processed_aggregated$Envelop))
    for (i in 1:EEE){
      inds = which(scan_processed_aggregated$Envelop==i)
      if (length(inds)>=2){ # At least 2 peaks in the envelop
        envelop = scan_processed_aggregated[inds,]
        results = annotate_envelop(envelop, transformation_list, IFL, ntheo = ntheo, baseline = baseline, min_overlap = min_overlap)
        
        scan_annotated = rbind.data.frame(scan_annotated,results$envelop.annotated)
        features_annotated = rbind.data.frame(features_annotated, results$features.annotated)
      }
    }} else {scan_annotated = features_annotated = NULL}
  
  if (!is.null(features_annotated)){
    if (nrow(features_annotated)>0){
      features_annotated = features_annotated[features_annotated$RESPONSE>0,]
      features_annotated1 = features_annotated[order(features_annotated$SCORE),]
      features_annotated = c()
      tmp_features = unique(features_annotated1$FEATURE)
      for (tf in tmp_features){
        valid = which(features_annotated1$FEATURE==tf)
        new_feature = features_annotated1[valid[1],]
        new_feature$RESPONSE = sum(new_feature$RESPONSE)
        features_annotated = rbind.data.frame(features_annotated, new_feature)
      }
    }
    
    valid = which(features_annotated$SCORE<=max_msigma & features_annotated$OC>= min_overlap & features_annotated$EXP_MMW_PPM<=max_mmw_ppm)
    features_annotated = features_annotated[valid,,drop=FALSE]
  }
  
  return(list(scan = scan_annotated, feature = features_annotated))
}

##############################################
####Annotate envelops by searching database###
##############################################

annotate_envelop<-function(envelop, ref_trans, IFL, ntheo, baseline, min_overlap){
  
  # Envelop annotation:
  
  envelop_annotated = c()
  
  sd1 = envelop
  tmp_cpd = tmp_formula = tmp_score = tmp_oc = tmp_cpd1 = rep(NA, nrow(sd1))
  tmp_feature = tmp_feature_envelop = tmp_feature_formula = c()
  tmp_feature_mmw = tmp_feature_amw = exp_feature_mmw = exp_feature_amw = c()
  exp_mmw_ppm = exp_amw_dev = exp_feature_ppm = exp_feature_dev = exp_feature_score = exp_feature_oc = exp_feature_response= c()
  
  coef1 = 1
  coef2 = 0
  
  tmp_scan = cbind.data.frame(sd1$MW, sd1$Response)
  avg_envelop = sum(tmp_scan[,1] * tmp_scan[,2]/sum(tmp_scan[,2]))
  max_envelop = max(tmp_scan[,1])
  min_envelop = min(tmp_scan[,1])

  valid = c()
  avg_dev = c()
  formula_valid = c()
  mSigma = NULL
  kt = c()
  for (k in 1:nrow(ref_trans)){
    dev_mass = abs(ref_trans$AVG.MW[k] - avg_envelop) 
    if (dev_mass<=5){
        RLF = useBRAIN(ListFormula1(ref_trans$FORMULA[k]), nrPeaks = ntheo)
        matched_nm = intersect(round(RLF$masses), round(tmp_scan[,1])) # Quickly check matched nominal mass
        if (length(matched_nm)>=ntheo*min_overlap){ # At least half of ntheo found
          kt = c(kt,length(matched_nm))
          valid = c(valid, k)
          avg_dev = c(avg_dev, dev_mass)
          formula_valid = c(formula_valid, ref_trans$FORMULA[k])
        }
   }}
  
  if (length(valid)>1){
    tnd = which(!duplicated(formula_valid))
    avg_dev = avg_dev[tnd]
    valid = valid[tnd]
  }

  # Only 1 possibility
  
  if (length(valid)==1){
    
    ifl = IFL[[valid]]
    theo_isotope <- useBRAIN(aC = ifl, nrPeaks = ntheo)
    mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope, ntheo, baseline)
    coef1 = 1
    coef2 = 0
  }
  
  # 2 possibilities:
  
  if (length(valid)>=2){

    tmp_r = rank(avg_dev) + rank(-kt)
    valid = valid[order(tmp_r)[1:2]]
    
    tmp_scan1 = tmp_scan
    tmp_scan1[,1] = tmp_scan1[,1] - 1.00726
    
    colnames(tmp_scan1) = c("mz", "intensity")
    ifl1 = IFL[[valid[1]]]
    ifl2 = IFL[[valid[2]]]
    theo_isotope1 <- useBRAIN(aC = ifl1, nrPeaks = ntheo)
    theo_isotope2 <- useBRAIN(aC = ifl2, nrPeaks = ntheo)
  
    optim <- try(deconvolution1(scan_df = tmp_scan1,
                        theor_ID_cmpd1 = theo_isotope1,
                        theor_ID_cmpd2 = theo_isotope2,
                        n_theor_peaks = ntheo,
                        expected_charge_range = 1,
                        matching_mass_accuracy = 0.15,
                        noise_threshold = 0,
                        deduplicate_fun = "max"
                  ), silent = F)
    
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
      mSigma = calcul_mix_mSigma(tmp_scan, theo_isotope1, theo_isotope2, coef1, coef2, ntheo, baseline)
    }
    if (coef1>0 & coef2==0){
      coef1 = 1
      valid = valid[1]
      mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope1, ntheo, baseline)
    }
    if (coef1==0 & coef2>0){
      coef2 = 1
      valid = valid[2]
      mSigma = calcul_imp_mSigma(tmp_scan, theo_isotope2, ntheo, baseline)
    }
  }
  
  # Output results:
  
  if (!is.null(mSigma)){
    if (!is.na(mSigma$score) & nrow(mSigma$exp_sp)>0){
     if (mSigma$score>0){
      if (!("score1" %in% names(mSigma)) & !is.null(mSigma$mono_mass)){
        coef1 = 1
        coef2 = 0
        inds = match(mSigma$exp_sp[,1], sd1$MW)
        inds = inds[!is.na(inds)]
        
        tmp_cpd[inds] = ref_trans$CPD[valid]
        tmp_formula[inds] = ref_trans$FORMULA[valid]
        tmp_score[inds] = round(mSigma$score,2)
        tmp_oc[inds] = round(mSigma$oc_score,2)
        tmp_cpd1[inds] = ref_trans$CPD[valid]
        
        tmp_feature = ref_trans$CPD[valid]
        tmp_feature_envelop = sd1$Envelop[1]
        tmp_feature_formula = ref_trans$FORMULA[valid]
        
        tmp_feature_mmw = round(mSigma$mono_mass_ref, 4) # Theoretical mono
        tmp_feature_amw = round(mSigma$avg_mass_ref, 4) # Theoretical average

        exp_feature_mmw = round(mSigma$mono_mass, 4)
        exp_feature_amw = round(mSigma$avg_mass, 4)
        exp_feature_ppm = round(mSigma$mono_ppm_dev,2)
        exp_feature_dev = round(mSigma$avg_mass_dev,2)
        exp_feature_score = round(mSigma$score,2)
        exp_feature_oc = round(mSigma$oc_score,2)
        
        exp_feature_response = round(sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp)]),0)
      } 
       
      if (("score1" %in% names(mSigma)) & !is.null(mSigma$mono_mass1)){

        inds1 = match(mSigma$exp_sp1[,1], sd1$MW)
        inds1 = inds1[!is.na(inds1)]
        tmp_cpd[inds1] = ref_trans$CPD[valid[1]]
        tmp_formula[inds1] = ref_trans$FORMULA[valid[1]]
        #tmp_score[inds1] = round(mSigma$score1,2)
        tmp_score[inds1] = round(mSigma$score,2)
  
        tmp_oc[inds1] = round(mSigma$oc_score1,2)
        
        inds2 = match(mSigma$exp_sp2[,1], sd1$MW)
        inds2 = inds2[!is.na(inds2)]
        tmp_cpd[inds2] = paste0(tmp_cpd[inds2], ":", ref_trans$CPD[valid[2]])
        tmp_formula[inds2] = paste0(tmp_formula[inds2], ":", ref_trans$FORMULA[valid[2]])
        #tmp_score[inds2] = paste0(tmp_score[inds2], ":", round(mSigma$score2,2))
        tmp_score[inds2] = paste0(tmp_score[inds2], ":", round(mSigma$score2,2))
        tmp_oc[inds2] =  paste0(tmp_oc[inds2], ":", round(mSigma$score,2))
        
        inds = unique(c(inds1, inds2))
        tmp_score[inds] = paste0(tmp_score[inds], "%", round(mSigma$score,2))
        tmp_oc[inds] = paste0(tmp_oc[inds], "%", round(mSigma$oc_score,2))

        tmp_cpd1[inds] =  paste0(ref_trans$CPD[valid[1]] , ":",ref_trans$CPD[valid[2]])
        
        tmp_feature = ref_trans$CPD[valid]
        tmp_feature_envelop = sd1$Envelop[1]
        tmp_feature_formula = ref_trans$FORMULA[valid]

        tmp_feature_mmw = round(c(mSigma$mono_mass_ref1, mSigma$mono_mass_ref2), 4) # Theoretical mono
        tmp_feature_amw = round(c(mSigma$avg_mass_ref1, mSigma$avg_mass_ref2),4) # Theoretical average
        
        exp_feature_mmw = round(c(mSigma$mono_mass1, mSigma$mono_mass2), 4)
        exp_feature_amw = round(c(mSigma$avg_mass1, mSigma$avg_mass2), 4)
        exp_feature_ppm = round(c(mSigma$mono_ppm_dev1, mSigma$mono_ppm_dev2),2)
        exp_feature_dev = round(c(mSigma$avg_mass_dev1, mSigma$avg_mass_dev2),2)
      #  exp_feature_score = round(c(mSigma$score1, mSigma$score2),2)
        exp_feature_score = round(c(mSigma$score, mSigma$score),2)
        exp_feature_oc = round(c(mSigma$oc_score1, mSigma$oc_score2),2)
        
        if (tmp_optim == "  there is no overlap"){
          tmp_r1 = sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp1)])
          tmp_r2 = sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp2)])
        } else {
          tmp_r1 = coef1*sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp)])
          tmp_r2 = coef2*sum(sd1$Response[which(sd1$MW %in% mSigma$exp_sp)])
        }
        exp_feature_response = round(c(tmp_r1, tmp_r2), 0)
      }
  }}}
  
  sd1$CPD = as.character(tmp_cpd)
  sd1$FORMULA = as.character(tmp_formula)
  sd1$SCORE = as.character(tmp_score)
  sd1$OC_SCORE  = as.character(tmp_oc)
  sd1$CPD = str_remove(sd1$CPD,"NA:")
  
  sd1$FORMULA = str_remove(sd1$FORMULA,"NA:")
  sd1$SCORE = str_remove(sd1$SCORE,"NA:")
  sd1$OC_SCORE = str_remove(sd1$OC_SCORE,"NA:")
  
  sd1$CPD1 = as.character(tmp_cpd1)
  sd1$COEF = paste0(round(coef1,2), ":", round(coef2,2))

  sd2 =  cbind.data.frame(FEATURE = tmp_feature, FORMULA = tmp_feature_formula, 
          THEO_MMW = tmp_feature_mmw, THEO_AMW = tmp_feature_amw,
          EXP_MMW = exp_feature_mmw, EXP_AMW = exp_feature_amw,
          EXP_MMW_PPM = exp_feature_ppm, EXP_AMW_DEV = exp_feature_dev, 
          SCORE = exp_feature_score, OC = exp_feature_oc,
          RESPONSE = exp_feature_response,
          Envelop = tmp_feature_envelop)
  
  return(list(envelop.annotated = sd1, features.annotated = sd2))
}

#####################################################
### Isotope evaluation functions single and mixed ###
#####################################################

calcul_imp_mSigma <- function(sp_deconvoluted, theo_isotope, ntheo, baseline){
  
  # Evaluate isotope patterns
  # sp_deconvoluted: m/z, intensity two column matrix
  
  # compute theoretical isotope distribution
  
  theo_deconvoluted = cbind.data.frame(theo_isotope$masses,theo_isotope$isoDistr)
  theo_deconvoluted[,2] = theo_deconvoluted[,2]/sum(theo_deconvoluted[,2])

  tmw = theo_deconvoluted[1,1] # Theoritical mono
  taw = sum(theo_deconvoluted[,1] * theo_deconvoluted[,2])
  
  sp_deconvoluted = as.data.frame(sp_deconvoluted)
  colnames(sp_deconvoluted) = colnames(theo_deconvoluted) = c("MW","I")
  NT = nrow(theo_deconvoluted)
  
  # Generate comparison vectors:

  em_list = rep(1, NT) # Dalton error from matched isotope
  res_list = rep(baseline/2, NT) # Response list from matched isotope
  exp_list = rep(0, NT) # Experimental mass, one-to-one matched to theoritical
  
  for (j in 1:NT){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted$MW[j])
    valid = which(errors<0.15)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list[j] = errors[valid]
      res_list[j] = sp_deconvoluted$I[valid]
      exp_list[j] = sp_deconvoluted$MW[valid]
    }
  }
  
  # Calculate chi_square and shape similarity score:
  
  tm_list = theo_deconvoluted$MW
  theo_list = theo_deconvoluted$I
  idx = which(res_list>baseline)
  
  if (length(idx)==0){return(NULL)}
  
 # oc1 = sum(res_list[idx])/sum(sp_deconvoluted[,2]) # Unnormalized
  #oc1 = length(idx)/nrow(sp_deconvoluted)
  #oc2 = sum(theo_deconvoluted$I[idx])/sum(theo_deconvoluted$I)
  oc2 = length(idx)/NT
  #oc_score = 2*(oc1*oc2)/(oc1+oc2)
  oc_score = oc2
    
  res_list = res_list/sum(res_list)
  theo_list = theo_list/sum(theo_list)
  
  chi_squa_score = sum(abs(res_list-theo_list)/theo_list)/NT
  
  # Compute evaluations:
    
  mono_mass = 0
  mono_ppm_dev = -1
  
  if (min(idx) == 1){ # Monoisotopic detected
    mono_mass = exp_list[1]
    mono_ppm_dev = abs(em_list[1])/tmw*1000000
  }
  
  avg_mass = sum(exp_list[idx] * res_list[idx]/sum(res_list[idx])) # Average molecular weight measured
  taw_bis = sum(tm_list * theo_list/sum(theo_list)) # Theoretical average in the filtered range
  avg_mass_dev = abs(avg_mass - taw_bis) # Average mass error in absolute!!!
  
  if (avg_mass>0){
    return(list(score = chi_squa_score, oc_score = oc_score, 
              mono_mass_ref= tmw, mono_mass =  mono_mass,
              avg_mass_ref = taw_bis, avg_mass = avg_mass, 
              mono_ppm_dev = mono_ppm_dev, avg_mass_dev = avg_mass_dev,
              exp_sp = cbind(Mass = exp_list, Res = res_list)))
  } 
}

calcul_mix_mSigma <- function(sp_deconvoluted, theo_isotope1, theo_isotope2, coef1, coef2, ntheo, baseline){
  
  # Evaluate isotope patterns
  # dat: m/z, intensity two column matrix
  
  # compute and weigh theoretical isotope distribution for two compounds
  
  #theo_isotope1 <- BRAIN::useBRAIN(aC = ifl1, nrPeaks = ntheo)
  #theo_isotope2 <- BRAIN::useBRAIN(aC = ifl2, nrPeaks = ntheo)
  
  theo_deconvoluted1 = cbind.data.frame(theo_isotope1$masses, theo_isotope1$isoDistr)
  theo_deconvoluted2 = cbind.data.frame(theo_isotope2$masses, theo_isotope2$isoDistr)
  theo_deconvoluted1[,2] = theo_deconvoluted1[,2]/sum(theo_deconvoluted1[,2])
  theo_deconvoluted2[,2] = theo_deconvoluted2[,2]/sum(theo_deconvoluted2[,2])
  colnames(sp_deconvoluted) = colnames(theo_deconvoluted1) = colnames(theo_deconvoluted2) = c("MW","I")
  
  tmw1 = theo_deconvoluted1[1,1] # Theoretical mono
  taw1 = sum(theo_deconvoluted1[,1] * theo_deconvoluted1[,2]) # Theoretical avg
  tmw2 = theo_deconvoluted2[1,1] 
  taw2 = sum(theo_deconvoluted2[,1] * theo_deconvoluted2[,2])
  
  # Mix two distributions, generate mixed theoritical distribution
  
  theo_deconvoluted1[,2] =  theo_deconvoluted1[,2]*coef1
  theo_deconvoluted2[,2] =  theo_deconvoluted2[,2]*coef2
  
  theo_mixed = rbind.data.frame(theo_deconvoluted1,theo_deconvoluted2)
  theo_mixed = theo_mixed[order(theo_mixed[,1]),]
  colnames(theo_mixed) = colnames(sp_deconvoluted) 
  theo_mixed$group = cut_mmw_list(theo_mixed$MW, theo_mixed$I, 0.9)$id
  theo_mass= aggregate(theo_mixed$MW, list(theo_mixed$group), FUN=mean) 
  theo_int= aggregate(theo_mixed$I, list(theo_mixed$group), FUN=sum) 
  theo_deconvoluted_mixed = cbind.data.frame(MW = theo_mass$x, I = theo_int$x)
  
  # Matching experimental spectra to theoretical:
  
  NT1 = nrow(theo_deconvoluted1)
  NT2 = nrow(theo_deconvoluted2)
  NT = nrow(theo_deconvoluted_mixed)
  
  em_list1 = rep(1, NT1) # Dalton error from matched isotope
  em_list2 = rep(1, NT2)
  em_list = rep(1, NT)
  res_list1 = exp_list1 = rep(baseline/2, NT1) # Response list from matched isotope
  res_list2 = exp_list2 = rep(baseline/2, NT2) # Experimental mass
  res_list = exp_list = rep(baseline/2, NT)
  
  for (j in 1:NT1){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted1$MW[j])
    valid = which(errors<0.15)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list1[j] = errors[valid]
      res_list1[j] = sp_deconvoluted$I[valid]
      exp_list1[j] = sp_deconvoluted$MW[valid]}
  }
  
  for (j in 1:NT2){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted2$MW[j])
    valid = which(errors<0.15)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list2[j] = errors[valid]
      res_list2[j] = sp_deconvoluted$I[valid]
      exp_list2[j] = sp_deconvoluted$MW[valid]}
  }
  
  for (j in 1:NT){
    errors = abs(sp_deconvoluted$MW - theo_deconvoluted_mixed$MW[j])
    valid = which(errors<0.15)
    if (length(valid)>1){valid = which.min(errors)}
    if (length(valid)==1){
      em_list[j] = errors[valid]
      res_list[j] = sp_deconvoluted$I[valid]
      exp_list[j] = sp_deconvoluted$MW[valid]}
  }
  
  # Calculate chi_square and shape similarity score:
  
  tm_list1 = theo_deconvoluted1$MW
  theo_list1 = theo_deconvoluted1$I
  tm_list2 = theo_deconvoluted2$MW
  theo_list2 = theo_deconvoluted2$I
  tm_list = theo_deconvoluted_mixed$MW
  theo_list = theo_deconvoluted_mixed$I
  idx1 = which(res_list1>baseline)
  idx2 = which(res_list2>baseline)
  idx = which(res_list>baseline)
  
 # oc1_1 = sum(res_list1[idx1])/sum(sp_deconvoluted[,2]) # Unnormalized
 # oc1_2 = sum(theo_list1[idx1])/sum(theo_list1)
  #oc1_1 = length(idx1)/nrow(sp_deconvoluted)
  oc1_2 =  length(idx1)/NT1
  #oc_score1 = 2*(oc1_1*oc1_2)/(oc1_1+oc1_2)
  oc_score1 = oc1_2
  
 # oc2_1 = sum(res_list2[idx2])/sum(sp_deconvoluted[,2]) # Unnormalized
 # oc2_2 = sum(theo_list2[idx2])/sum(theo_list2)
  #oc2_1 = length(idx2)/nrow(sp_deconvoluted)
  oc2_2 =  length(idx2)/NT2
 # oc_score2 = 2*(oc2_1*oc2_2)/(oc2_1+oc2_2)
  oc_score2 = oc2_2
  
  #oc_1 = sum(res_list[idx])/sum(sp_deconvoluted[,2]) # Unnormalized
  #oc_2 = sum(theo_list[idx])/sum(theo_list)
  #oc_1 = length(idx)/nrow(sp_deconvoluted)
  oc_2 =  length(idx)/NT
  #oc_score = 2*(oc_1*oc_2)/(oc_1+oc_2)
  oc_score = oc_2
  
  res_list1 = res_list1/sum(res_list1)
  res_list2 = res_list2/sum(res_list2)
  res_list = res_list/sum(res_list)
  
  theo_list1 = theo_list1/sum(theo_list1)
  theo_list2 = theo_list2/sum(theo_list2)
  theo_list = theo_list/sum(theo_list)
  
  chi_squa_score1 = sum(abs(res_list1-theo_list1)/theo_list1)/NT1
  chi_squa_score2 = sum(abs(res_list2-theo_list2)/theo_list2)/NT2
  chi_squa_score = sum(abs(res_list-theo_list)/theo_list)/NT
  
  # Compute deviations:
  
  mono_mass1 = 0
  mono_ppm_dev1 = -1
  if (min(idx1)==1){
     # Monoisotopic detected
    mono_mass1 = exp_list1[1]
    mono_ppm_dev1 = abs(em_list1[1])/tmw1*1000000
  }
  
  mono_mass2 = 0
  mono_ppm_dev2 = -1
  if (min(idx2)==1){ # Monoisotopic detected
      mono_mass2 = exp_list2[1]
      mono_ppm_dev2 = abs(em_list2[1])/tmw2*1000000
  }

  #avg_mass1 = sum(exp_list1[idx1]*res_list1[idx1]/sum(res_list1[idx1])) # Average molecular weight measured
  #avg_mass2 = sum(exp_list2[idx2]*res_list2[idx2]/sum(res_list2[idx2]))
  #avg_mass = sum(exp_list[idx]*res_list[idx]/sum(res_list[idx])) 
  
  avg_mass1 = sum(exp_list1[idx1]*theo_list1[idx1]/sum(theo_list1[idx1])) # Average molecular weight measured
  avg_mass2 = sum(exp_list2[idx2]*theo_list2[idx2]/sum(theo_list2[idx2]))
  avg_mass = sum(exp_list[idx]*theo_list[idx]/sum(theo_list[idx])) 
  
  taw_bis1 = sum(tm_list1*theo_list1/sum(theo_list1)) # Theoritical average in the filtered range
  taw_bis2 = sum(tm_list2*theo_list2/sum(theo_list2)) 
  taw_bis = sum(tm_list*theo_list/sum(theo_list)) 
  
  avg_mass_dev1 = abs(avg_mass1 - taw_bis1) # Average mass error in absolute!!!
  avg_mass_dev2 = abs(avg_mass2 - taw_bis2) 
  avg_mass_dev = abs(avg_mass - taw_bis) 
  
  if (avg_mass1>0 & avg_mass2>0){
    return(list(score1 = chi_squa_score1, score2 = chi_squa_score2, score = chi_squa_score,
        oc_score1 = oc_score1, oc_score2 = oc_score2, oc_score = oc_score,
        mono_mass_ref1 = tmw1, avg_mass_ref1 = taw_bis1, mono_mass1 =  mono_mass1, avg_mass1 = avg_mass1,
        mono_mass_ref2 = tmw2, avg_mass_ref2 = taw_bis2, mono_mass2 =  mono_mass2, avg_mass2 = avg_mass2,
        avg_mass_ref = taw_bis, avg_mass = avg_mass,  avg_mass_dev = avg_mass_dev,
        mono_ppm_dev1 = mono_ppm_dev1, avg_mass_dev1 = avg_mass_dev1,
        mono_ppm_dev2 = mono_ppm_dev2, avg_mass_dev2 = avg_mass_dev2,
        exp_sp1 = cbind(Mass = exp_list1, Res = res_list1),
        exp_sp2 = cbind(Mass = exp_list2, Res = res_list2),
        exp_sp = cbind(Mass = exp_list, Res = res_list)))
  }
}

#################################
### Prepare impurity database ###
#################################

expand_transformation_list<-function(formula_flp, transformation_list, mdb){
  
  #transformation_list = read.csv("Transformation_list_jennifer.txt", sep = "\t", stringsAsFactors = F)
  # Oxaluric Acid (-20), Unlock nucleic acid, Deamination (+1Da), 2' Fluorination (+19)
  # Spiroiminodihydantion (+32), (S)-constrained ethyl (+26), 2'-Methoxy (O-methyl)(+31)
  # 2'-O-ally(40), Acetylation of Cytosine (+42), Locked nucleic acid (44),ADMF/GDMF (57), 2' Methoxyethyl(58)
  
  if (is.null(mdb)){
    IFL = calcul_imp_formula(formula_flp, transformation_list)
    amwFLP = calculateAverageMass(ListFormula1(formula_flp))
    mmwFLP = calculateMonoisotopicMass(ListFormula1(formula_flp))
  } else {
    transformation_list = NULL
    IFL = calcul_imp_formula("N/A", mdb)
    amwFLP = mmwFLP = 0
    transformation_list$CPD = mdb$CPD
  } 
  
  NIF = length(IFL)
  
  if (length(NIF)>0){
    scan_annotated = c()
    features_annotated = c()
  
    for (i in 1:NIF){
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
    transformation_list = data.frame(transformation_list)
  }
  return(list(IFL = IFL, transformation_list = transformation_list))
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


