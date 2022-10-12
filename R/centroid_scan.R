#' Centroiding mass spectra
#'
#' The function creates centroid-mode mass spectrum if profile mode mass spectrum is detected
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom MSnbase pickPeaks
#' @export

centroid_scan<-function(scan0){
  
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
