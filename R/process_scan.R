#' Converting multi-charged oligonucleotide mass spectra to true molecular weight scale
#'
#' The function determines charges state based on peak spacing. It also creates a deconvoluted molecular weight spectra by combining multiple charge state of the same isotopic species.
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#'
#' @export
#'
##################################
######## Main Function ###########
##################################

process_scan<-function(test.scan, polarity = c("Positive", "Negative"),  MSMS = F, baseline = 100,
                       min_charge = 3, max_charge = 12, min_mz = 500, max_mz = 1500, min_mw = 4000, max_mw = 12000, 
                       mz_error = 0.02, mw_gap = 5){
  
  options(stringsAsFactors = F)
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


###############################
# Find high charge state peaks#
###############################

process_scan_high_charge_bis<-function(scan1, ref_charge, mz_error){
  
  colnames(scan1) = c("Mass", "Response")
  scan1 = data.matrix(scan1)
  NS = nrow(scan1)
  NC = length(ref_charge)
  ref_dis = 1/ref_charge # Difference between monoisotopic peaks for different charge state
  
  exp_mz = scan1[,1]
  
  # Validate charge state by looking at 0.5 Da before/after:
  
  scan2 = c()
  
  for (j in 1:NS){
    
    lookup_range = which(exp_mz>=exp_mz[j]-0.5 & exp_mz<=exp_mz[j]+0.5)
    lookup_range = max(1, j-3): min(j+3, NS)
    mzl = exp_mz[lookup_range]
    dist_mzl = data.matrix(dist(mzl))
    dist_mzl = dist_mzl[lower.tri(dist_mzl, diag = FALSE)]
    dev_charge = sapply(1:NC, function(x) abs(dist_mzl - ref_dis[x]))
    tmp_matched = apply(dev_charge, 1, function(x) {ifelse(sum(x<=mz_error)>0, 1, 0)})
    charge_mz = ref_charge[as.matrix(apply(dev_charge, 1, which.min))]
    charge_mz = charge_mz*tmp_matched
    charge_mz = charge_mz[charge_mz>0]
    
    if (length(charge_mz)>0){
      charge_count = sort(table(charge_mz), decreasing =TRUE)
      if (charge_count[1]>2){
        best_charge = as.numeric(names(charge_count[1]))
        temp_dat = scan1[j,,drop=FALSE]
        temp_dat =  cbind.data.frame(temp_dat, z = best_charge)
        scan2 = rbind(scan2, temp_dat)
      }
    }
    
    # Output:
    
    if (!is.null(scan2)){
      colnames(scan2) = c("Mass", "Response", "z")
      scan2 = data.frame(scan2)
    }
  }
  return(scan2)
}

##############################
# Find low charge state peaks#
##############################

process_scan_low_charge<-function(scan1, ref_charge, mz_error, baseline){
  
  # The function must be used after high mass determination
  
  scan3 = c()
  
  # Determine double charge
  
  if (2 %in% ref_charge & nrow(scan1)>3){
    
    exp_mz = scan1[,1]
    dist_mz = data.matrix(dist(exp_mz))
    
    tmp_inds = which(abs(dist_mz - 1/2)<=mz_error, arr.ind = T) # Return pairs of peaks
    all_pairs = tmp_inds[tmp_inds[,1] - tmp_inds[,2]<0,,drop=FALSE]
    NP = nrow(all_pairs)
    
    if (NP>0){
      
      for (j in 1:NP){
        
        tmp_range = all_pairs[j,1]:all_pairs[j,2]
        
        Segment = scan1[tmp_range,,drop=FALSE]
        SGR = Segment$Response[1]/Segment$Response[nrow(Segment)]
        
        # Check shape
        
        if (SGR>1.1 & nrow(Segment)<6){
          Segment = cbind.data.frame(Segment, z = 2)
          scan3 = rbind.data.frame(scan3, Segment)
        }
      }
      scan1 = scan1[!(scan1$Mass %in% scan3$Mass),,drop = FALSE]
    }
  }
  
  # Determine mono charge
  
  if (1 %in% ref_charge & nrow(scan1)>3){
    
    exp_mz = scan1[,1]
    dist_mz = data.matrix(dist(exp_mz))
    
    tmp_inds = which(abs(dist_mz - 1)<=mz_error, arr.ind = T) # Return pairs of peaks
    all_pairs = tmp_inds[tmp_inds[,1] - tmp_inds[,2]<0,,drop=FALSE]
    NP = nrow(all_pairs)
    
    if (NP>0){
      
      for (j in 1:NP){
        
        tmp_range = all_pairs[j,1]:all_pairs[j,2]
        
        Segment = scan1[tmp_range,,drop=FALSE]
        SGR = Segment$Response[1]/Segment$Response[nrow(Segment)]
        
        # Check shape
        
        if (SGR>1.1 & nrow(Segment)<6){
          Segment = cbind.data.frame(Segment, z = 1)
          scan3 = rbind.data.frame(scan3, Segment)
        }
      }}
    
    scan3 = scan3[!duplicated(scan3$Mass),,drop=FALSE]
    scan1 = scan1[!(scan1$Mass %in% scan3$Mass),,drop = FALSE]
  }
  
  # Determine unknown charge
  
  if (1 %in% ref_charge & nrow(scan1)>1){
    Segment = scan1[scan1$Response>baseline*10,,FALSE]
    if (nrow(Segment)>0){
      Segment = cbind.data.frame(Segment, z = 1)
      scan3 = rbind.data.frame(scan3, Segment)
      scan3 = scan3[!duplicated(scan3$Mass),,drop=FALSE]
    }}
  
  return(scan3)
}

##################################################
# Aggregate features with different charge states#
##################################################

process_aggregation<-function(scan_processed){
  
  # The function aggregate based on temporary features:
  
  tmp_feature = scan_processed$tmp_feature
  all_features = unique(tmp_feature)
  NF = length(all_features)
  new_mass = rep("0", NF)
  new_Res = rep(0, NF)
  new_z = rep("0", NF)
  new_MW = rep(0, NF)
  
  for (i in 1:NF){
    
    valid = which(tmp_feature == all_features[i])
    tmp_scan = scan_processed[valid,,drop = FALSE]
    
    masslist = sort(as.numeric(unlist(sapply(as.numeric(tmp_scan$Mass), function(x) strsplit(as.character(x), ":")[[1]]))))
    new_mass[i] = paste0(round(masslist,4), collapse = ":")
    new_Res[i] = sum(tmp_scan$Response)
    
    zlist = sort(unique(as.numeric(unlist(sapply(as.numeric(tmp_scan$z), function(x) strsplit(as.character(x), ":")[[1]])))))
    new_z[i] = paste0(zlist, collapse = ":")
    new_MW[i] = round(mean(tmp_scan$MW),3)
  }
  
  aggregated_scan = cbind.data.frame(Mass = new_mass, Response = new_Res, z = new_z, MW = new_MW)
  
  return(aggregated_scan)
}

#######################
### Other functions####
#######################

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

cut_mmw_list1<-function(mwlist, intlist, mw_window){
  
  # Used to cluster/separate envelops!  
  
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
    
    if ((mwlist[k] - min_mw > mw_window & mwlist[k] - best_mw > mw_window/2) || (mwlist[k] - max_mw >=1.2)){
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
