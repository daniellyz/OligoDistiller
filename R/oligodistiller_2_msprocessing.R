library(BRAIN)
library(OrgMassSpecR)
options(stringsAsFactors = F)

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
  
  # Validate charge state by looking at 3 neighboring peak before/after:
  
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
### Helper functions###
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
