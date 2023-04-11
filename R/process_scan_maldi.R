#' Converting maldi oligonucleotide mass spectra
#'
#' The function processes and envelops maldi scan
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#'
#' @export
#'
##################################
######## Main Function ###########
##################################

process_scan_maldi<-function(test.scan, polarity = c("Positive", "Negative"),  baseline = 1000,
                       min_mz = 100, max_mz = 2500,
                       mw_gap = 1.1, mw_window = 10){
  
  options(stringsAsFactors = F)
  
  scan0 = test.scan
  if (is.null(scan0)){scan0 = matrix(0,2,2)}
  if (nrow(scan0)<=5){scan0 = matrix(0,2,2)}
  
  scan_processed_aggregated = c()
  scan_processed = c()  
  
  # Pre-processing the scan:
  
  if (nrow(scan0)>5){
    scan0 = scan0[order(scan0[,1]),]
    scan0 = cbind(scan0[,1], scan0[,2])
    scan0 = apply(scan0, 2, as.numeric)
    scan0 = data.frame(scan0)
    colnames(scan0) = c("Mass", "Response")
    scan_processed = scan0[scan0$Response>baseline & scan0$Mass>=min_mz & scan0$Mass<=max_mz,,drop =FALSE]
  }
    
  if (is.null(scan_processed)){scan_processed = matrix(0,2,2)}
    
  if (nrow(scan_processed)>5){
      
      if (polarity == "Positive"){
      
        scan_processed = cbind.data.frame(Mass = scan_processed[,1], Response = scan_processed[,2], z=1, MW = scan_processed[,1]-1.00726)
      }
    
     if (polarity == "Negative"){
      
       scan_processed = cbind.data.frame(Mass = scan_processed[,1], Response = scan_processed[,2], z=1, MW = scan_processed[,1]+1.00726)
     }
    
        
      mmw_cutted = cut_mmw_list1(scan_processed$MW, scan_processed$Response, mw_gap, mw_window)
      
      scan_processed_aggregated = scan_processed
      scan_processed_aggregated$Envelop = mmw_cutted$id
    
    }
    return(list(scan_processed = scan_processed, scan_processed_aggregated = scan_processed_aggregated))
}

    
#######################
### Other functions####
#######################
    
cut_mmw_list1<-function(mwlist, intlist, mw_gap, mw_window){
      
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
        ratio_int = intlist[k]/intlist[k-1]
        if (k>2){ratio_int0 =  intlist[k-1]/intlist[k-2]} else {ratio_int0 = 2}
        
        if ((mwlist[k] - best_mw > mw_window/2*1.1  & ratio_int>1.1 & ratio_int0>1.1) 
            || (mwlist[k] - min_mw > mw_window*1.1  & ratio_int>1.1 & ratio_int0>1.1) 
            || (mwlist[k] - max_mw >mw_gap*1.1)){
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
    