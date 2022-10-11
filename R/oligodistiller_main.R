
library(MSnbase)
library(tools)
library(pracma)
library(plyr)
library(stringr)
source("oligodistiller_2_msprocessing.R")
source("oligodistiller_3_annotation.R")
source("oligodistiller_3_unmixing.R")

#################
## Centroiding ##
#################

centroid_scan<-function(scan0){
  
  # Detect profile mode data
  
  centroided = T
  scan0 = scan0[which(!is.na(scan0[,1])),]
  scan0 = scan0[which(!is.na(scan0[,2])),]
  scan0[which(scan0[,2]<0),2] = 0
  
  tmp_mz = scan0[,1]
  tmp_mz1 = c(0, tmp_mz)
  tmp_mz2 = c(tmp_mz, 0)
  tmp_md = tmp_mz2 - tmp_mz1
  if (sum(abs(tmp_md)<0.9, na.rm = T)>=5){centroided = F}
  
  if (!centroided){
    tmp_sp <- new("Spectrum1",
                  intensity = scan0[,2],
                  mz = scan0[,1],
                  centroided = FALSE)
    scan1 <- pickPeaks(tmp_sp)
    scan1 = cbind(scan1@mz, scan1@intensity)
  } else {scan1 = scan0}
  return(scan1)
}

############################################################
## Determine charge, molecular weight and isotope envelop###
############################################################

process_scan<-function(test.scan, polarity = c("Positive", "Negative"),  MSMS = F, baseline = 100,
          min_charge = 3, max_charge = 12, min_mz = 500, max_mz = 1500, min_mw = 4000, max_mw = 12000, 
          mz_error = 0.02, mw_gap = 5){
  
  if (min_charge>max_charge){min_charge = max_charge-1}
  
  ref_charge_high = c()
  ref_charge_low = c()
  if (max_charge>=3){ref_charge_high = max(3, min_charge):max_charge}
  if (min_charge<3){ref_charge_low = min_charge:2}
  
  scan0 = test.scan
  if (is.null(scan0)){scan0 = matrix(0,2,2)}
  if (nrow(scan0)<=5){scan0 = matrix(0,2,2)}
  
  scan_processed_aggregated = c()
  scan_processed = c()
  
  # Pre-processing the scan:
  
  if (nrow(scan0)>5){
    scan0 = cbind(scan0[,1], scan0[,2])
    scan0 = apply(scan0, 2, as.numeric)
    scan0 = data.frame(scan0)
    colnames(scan0) = c("Mass", "Response")
    scan0 = scan0[scan0$Response>baseline & scan0$Mass>=min_mz & scan0$Mass<=max_mz,,drop =FALSE]
  }
  
  # High charge state processing:
  
  if (nrow(scan0)>5 & length(ref_charge_high)>0){
    scan1 = process_scan_high_charge_bis(scan0, ref_charge_high, mz_error)
    scan_processed = rbind.data.frame(scan_processed, scan1)
    scan0 = scan0[!(scan0$Mass %in% scan_processed$Mass),,drop = FALSE] # Filter original scan
  }
  
  # Low charge state processing:
  
  if (nrow(scan0)>5 & length(ref_charge_low)>0){
    scan1 = process_scan_low_charge(scan0, ref_charge_low, mz_error, baseline)
    scan_processed = rbind.data.frame(scan_processed, scan1)
  }
  
  # Construct and filter scan:
  # Minimum 2 isotopic species
  
  if (is.null(scan_processed)){scan_processed = matrix(0,2,2)}
  
  if (nrow(scan_processed)>5){
    
    if (polarity == "Positive"){scan_processed$MW = scan_processed$Mass*scan_processed$z - scan_processed$z*1.00726}
    if (polarity == "Negative"){scan_processed$MW = scan_processed$Mass*scan_processed$z + scan_processed$z*1.00726}
    scan_processed = scan_processed[order(scan_processed$MW),,drop=FALSE]
    scan_processed = scan_processed[which(scan_processed$MW>=min_mw & scan_processed$MW<=max_mw),,drop=FALSE]
    
    if (nrow(scan_processed)>5){
      
      mmw_cutted = cut_mmw_list(scan_processed$MW, scan_processed$Response, 1)
      scan_processed$tmp_feature = mmw_cutted$id
      scan_processed$MW = mmw_cutted$mw
      
      if (MSMS){frequent.feature = as.numeric(which(table(scan_processed$tmp_feature)>=1))}
      if (!MSMS){frequent.feature = as.numeric(which(table(scan_processed$tmp_feature)>1))}
      scan_processed = scan_processed[which(scan_processed$tmp_feature %in% frequent.feature),,drop=FALSE]
      
      if (nrow(scan_processed)>3){
        
        scan_processed_aggregated = process_aggregation(scan_processed)
        scan_processed_aggregated$Envelop = cut_mmw_list1(scan_processed_aggregated$MW, scan_processed_aggregated$Response, mw_gap)$id
      }}
  }
  
  return(scan_processed_aggregated)
}

################################
###Deisotoping & Annotation ####
################################

annotate_scan_targeted<-function(scan_processed_aggregated, formula_flp = "C192H239O117N73P18S4F8", transformation_list = NULL, ntheo = 12){
  
  # Calculate impurity formula:
  
  flp_mw_avg = MolecularWeight(formula = ListFormula1(formula_flp))
  
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
  
  # Annotate envelop
  
  if (!is.null(scan_processed_aggregated)){
    EEE = max(as.numeric(scan_processed_aggregated$Envelop))
    for (i in 1:EEE){
      inds = which(scan_processed_aggregated$Envelop==i)
      if (length(inds)>=3){ # At least 3 peaks in the envelop
        envelop = scan_processed_aggregated[inds,]
        results = annotate_envelop(envelop, transformation_list, IFL, ntheo = 12)
        scan_annotated = rbind.data.frame(scan_annotated,results$envelop.annotated)
        features_annotated = rbind.data.frame(features_annotated, results$features.annotated)
      }
    }} else {scan_annotated = features_annotated = NULL}
  
  if (!is.null(features_annotated)){
    features_annotated$MMW = round(as.numeric(features_annotated$AMW), 4)
    features_annotated$SCORE = round(as.numeric(features_annotated$SCORE), 2)
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
  return(list(scan = scan_annotated, feature = features_annotated))
}

annotate_scan_targeted_wo_flp<-function(scan_processed_aggregated, mdb = NULL, ntheo = 12){
  
  # Calculate impurity formula:
  
  IFL = calcul_imp_formula("N/A", mdb)

  NIF = length(IFL)
  scan_annotated = c()
  features_annotated = c()
  mdb$AVG.MW = sapply(IFL, BRAIN::calculateAverageMass)
  mdb$MONO.MW = sapply(IFL, BRAIN::calculateMonoisotopicMass)
  
  # Annotate envelop
  
  if (!is.null(scan_processed_aggregated)){
    EEE = max(as.numeric(scan_processed_aggregated$Envelop))
    for (i in 1:EEE){
      inds = which(scan_processed_aggregated$Envelop==i)
      if (length(inds)>=3){ # At least 3 peaks in the envelop
        envelop = scan_processed_aggregated[inds,]
        results = annotate_envelop(envelop, mdb, IFL, ntheo = 12)
        scan_annotated = rbind.data.frame(scan_annotated,results$envelop.annotated)
        features_annotated = rbind.data.frame(features_annotated, results$features.annotated)
      }
    }} else {scan_annotated = features_annotated = NULL}
  
  if (!is.null(features_annotated)){
    features_annotated$AMW = round(as.numeric(features_annotated$AMW), 4)
    features_annotated$SCORE = round(as.numeric(features_annotated$SCORE), 2)
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
  
  return(list(scan = scan_annotated, feature = features_annotated))
}

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

####################################
###Reconstruct original spectra ####
####################################

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