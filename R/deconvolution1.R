#' Deconvoluting a mixed envelop to two pure components
#'
#' The function decomposes a mixed envelop into two theoritical isotope envelopes. It quantifies the two components by decomposing.
#'
#' @param scan_df Data frame with two columns, mz and intensity, representing all peaks in a single mass scan at one selected retention time
#' @param ac_cmpd1 List used by the BRAIN algorithm to compute theoretical isotope distributions. The lists have to follow this format: list(C = 1, H = 2, F = 3, N = 4, O = 5, P = 6, S = 7)), where numbers 1-7 should be replaced with actual values
#' @param ac_cmpd2 Same as previous input, another expected theoritical envelop list
#' @param n_theor_peaks Numeric. How many theoretical peaks should be returned by BRAIN
#' @param expected_charge_range Numeric, which charge state should be analysed? e.g. 5:11
#' @param matching_mass_accuracy Numeric, relative (in ppm) mass tolerance for assigning observed peaks to theoretical peaks
#' @param noise_threshold Numeric, keep for the analysis only those peaks with intensities greater or equal than this threshold
#' @param extraction_mass_accuracy Numeric, 50 by default. mass tolerance (in ppm) for isotopic peak cluster extraction from the entire mass spectrum
#' @param deduplicate_fun Character, "max" by default. How to deduplicate multiple observed peaks assigned to one theoretical peak - "max" or "sum"
#' 
#' @return
#' \itemize{
#'    \item{by_charge:}{ list, with elements like z5, z6, ..., each storing estimate (estimated proportion of cmpd2, 0-1 scale),
#'                    se (standard error of the proportion estimate), p_value (null hypothesis: estimate=0), 
#'                    significance (1 = significant, 0 = insignificant), mpcse (pearson chi-square goodness-of-fit, similar to the pointless4dna paper)}
#'    \item{by_charge:}{ not yet available}
#'  }
#'  
#' @author Piotr Prostko, \email{piotr.prostko@uhasselt.be}
#' 
#' @importFrom tibble tibble
#' @importFrom purrr map_dfc pmap set_names
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_c
#' @importFrom dplyr bind_cols bind_rows mutate select summarise group_by filter ungroup lag
#' 
#' @export
#'

deconvolution1 <- function(scan_df, theor_ID_cmpd1, theor_ID_cmpd2, n_theor_peaks, expected_charge_range, matching_mass_accuracy, noise_threshold, extraction_mass_accuracy = 50, deduplicate_fun = "max", proton_mass = 1.007276466621) {

  # Without giving mass value  
  # remove noisy peaks
  
  scan_df <- scan_df[scan_df[,2] >= noise_threshold,]
  
  # Da difference between template spectra
  
  template_mass_diff <- round(theor_ID_cmpd1$monoisotopicMass-theor_ID_cmpd2$monoisotopicMass)
  
  # extract isotopic envelopes from the observed mass spectrum
  
  isotopes <- extract_isotopes(scan_df, theor_ID_cmpd1, theor_ID_cmpd2, proton_mass, expected_charge_range)
  charge_df <- isotopes$charge_df
  isotopes <- isotopes$data
  
  # compute mass to charge ratios for theoretical IDs
  theor_ID_cmpd1 <- set_names(charge_df$z, str_c(charge_df$z)) %>% 
    map_dfc(~(theor_ID_cmpd1$masses-.x*proton_mass)/.x) %>% 
    bind_cols(intensity = theor_ID_cmpd1$isoDistr) %>% 
    pivot_longer(cols = -"intensity", names_to = "z", values_to = "mz") %>% 
    # mutate(spectrum = "cmpd1", z = as.numeric(z), .before = everything())
    mutate(spectrum = "cmpd1", z = as.numeric(z)) %>% 
    select(spectrum, intensity, z, mz) # Piotr's update: this order of columns is required
  
  theor_ID_cmpd2 <- set_names(charge_df$z, str_c(charge_df$z)) %>% 
    map_dfc(~(theor_ID_cmpd2$masses-.x*proton_mass)/.x) %>% 
    bind_cols(intensity = theor_ID_cmpd2$isoDistr) %>% 
    pivot_longer(cols = -"intensity", names_to = "z", values_to = "mz") %>% 
    # mutate(spectrum = "cmpd2", z = as.numeric(z), .before = everything())
    mutate(spectrum = "cmpd2", z = as.numeric(z)) %>% 
    select(spectrum, intensity, z, mz) # Piotr's update: this order of columns is required
  
  # the spectrum shifted to the left should be used for matching the observed peaks with theoretical peaks
  if (template_mass_diff > 0) { # cmpd1 mono mass - cmpd2 mono mass
    ID_match <- theor_ID_cmpd2
  } else ID_match <- theor_ID_cmpd1
  
  mz_tmp <- ID_match$mz
  alignment_mz_bounds <- data.frame(z = ID_match$z, mz_lb = mz_tmp - mz_tmp * matching_mass_accuracy / 10^6, mz_ub = mz_tmp + mz_tmp * matching_mass_accuracy / 10^6)
  
  # split data by charge
  isotopes_z_split <- split(isotopes, isotopes$z)  
  theor_ID_cmpd1_z_split <- split(theor_ID_cmpd1, theor_ID_cmpd1$z)
  theor_ID_cmpd2_z_split <- split(theor_ID_cmpd2, theor_ID_cmpd2$z)
  alignment_mz_bounds_z_split <- split(alignment_mz_bounds, alignment_mz_bounds$z)
  
  results <- vector("list", 2) %>% set_names("by_charge", "combined_charge")
  results$by_charge <- vector("list", length(isotopes_z_split)) %>% set_names(paste0("z", names(isotopes_z_split)))
  
  for (z_ind in seq_along(isotopes_z_split)){
    # match (in the mz dimension) observed peaks to the theoretical isotope distribution (given by ID_match)
    data <- align_peaks(data = isotopes_z_split[[z_ind]][,-1], 
                        bound = alignment_mz_bounds_z_split[[z_ind]][,-1])
    # if multiple observed peaks correspond to one theoretical peaks, take the most abundant one
    data <- deduplicate_peaks(data, deduplicate_fun)
    TIC <- sum(data$intensity, na.rm = TRUE)
    # sum=1 normalization
    data$intensity <- data$intensity/TIC
    # create design matrix (with leading and trailing zeroes that correspond to the mass difference between templates)
    design_matrix <- get_design_matrix(X = cbind(cmpd1 = theor_ID_cmpd1_z_split[[z_ind]]$intensity, 
                                                 cmpd2 = theor_ID_cmpd2_z_split[[z_ind]]$intensity),
                                       match_index = data$theor_ind, mass_diff = template_mass_diff, normalize = TRUE)
    mod <- nls(I(mix - cmpd1) ~ theta * I(cmpd2 - cmpd1), start = c("theta" = 0.5), algorithm = "port", lower = 0, upper = 1,
               data = cbind(mix = data$intensity, design_matrix))
    mod_sum <- summary(mod)
    
    mod_fitted <- fitted(mod) + design_matrix$cmpd1
    N <- TIC / sum(mod_fitted)
    tmp_ind <- is.finite(data$intensity) & is.finite(mod_fitted)
    
    tmp_mpcse <- TIC*(data$intensity[tmp_ind] - mod_fitted[tmp_ind])^2/mod_fitted
    mpcse <- mean(tmp_mpcse[is.finite(tmp_mpcse)], na.rm = TRUE)
    
    results$by_charge[[z_ind]] <- c(estimate = mod_sum$coefficients[, 1], mpcse = mpcse)
  }
  return(results)
}

# example -----------------------------------------------------------------

#source("exmapledata4Youzhong.R")
#optim <- deconvolution(scan_df = dat_dev,
#             ac_cmpd1 = list(C = 210, H = 262, F = 9, N = 75, O = 132, P = 20, S = 4),
#             ac_cmpd2 = list(C = 210, H = 263, F = 8, N = 75, O = 133, P = 20, S = 4),
#             n_theor_peaks = 30,
#             expected_charge_range = 5:11,
#             matching_mass_accuracy = 10,
#             noise_threshold = 1000)


# bm <- bench::mark(deconvolution(scan_df = dat_dev,
#                                 ac_cmpd1 = list(C = 210, H = 262, F = 9, N = 75, O = 132, P = 20, S = 4),
#                                 ac_cmpd2 = list(C = 210, H = 263, F = 8, N = 75, O = 133, P = 20, S = 4),
#                                 n_theor_peaks = 30,
#                                 expected_charge_range = 5:11,
#                                 matching_mass_accuracy = 10,
#                                 noise_threshold = 1000),
#             iterations = 100)
# bm


######################
### Helper functions##
######################

norm <- function(x) x/sum(x, na.rm = TRUE)

align_peaks <- function(data, bound){
  # match the observed mixture and theoretical spectra over the m/z dimension
  # do this using only one of the theoretical spectra
  mz <- data$mz  
  ind <- map(mz, ~ (.x >= bound[,1] & .x <= bound[,2]))
  ind2 <- ind %>% map(~ any(.x)) %>% unlist()
  ind3 <- ind %>% map(~which(.x)) %>% map(~ifelse(length(.x)==0, 0, .x)) %>% unlist()
  return(data.frame(data[ind2, ], theor_ind=ind3[ind2]))
}

deduplicate_peaks <- function(data, method){
  method.arg <- match.arg(method, c("max", "sum"))
  if (method.arg == "max"){
    # simple peak picking - take the max if multiple peaks correspond to one theoretical peak
    data <- data %>% 
      group_by(theor_ind) %>%
      filter(intensity == max(intensity)) %>% 
      ungroup() 
  } else if (method.arg == "sum"){
    data <- data %>% 
      group_by(theor_ind) %>%
      summarise(intensity = sum(intensity)) %>% 
      ungroup() 
  }
  return(data)
}

get_design_matrix <- function(X, match_index, mass_diff, normalize = TRUE) {
  if (abs(mass_diff) != floor(abs(mass_diff))) stop("mass_diff is not integer")
  if (abs(mass_diff) >= length(match_index)) stop("there is no overlap")
  if (length(match_index) < 2) stop("The length of match_index should be >= 2")
  
  if (mass_diff > 0) {
    cmpd2 = c(X[1 : (nrow(X)-mass_diff), "cmpd2"], rep(0, mass_diff))
    cmpd2 = cmpd2[match_index]
    cmpd1 = dplyr::lag(X[, "cmpd1"], n = mass_diff, default = 0)[match_index]
  } else {
    cmpd2 = dplyr::lag(X[, "cmpd2"], n = abs(mass_diff), default = 0)[match_index]
    cmpd1 = c(X[1 : (length(match_index)-abs(mass_diff)), "cmpd1"], rep(0, abs(mass_diff)))
    cmpd1 = cmpd1[match_index]
  }
  out <- tibble(cmpd1 = cmpd1, cmpd2 = cmpd2)
  if (normalize) out <- map_dfc(out, ~norm(.x))
  return(out)
  
}

# extract isotopic envelopes from the observed mass spectrum
# the algorithm for the peak extractions goes as follows:
# 1) for the two theoretical ID, determine the number of peaks that sum up to the probability_coverage value (e.g. 0.99)
# 2) combine mass values corresponding to the selected peaks in the step above and compute their range
# 3) this mass range is used in computing lower and upper mz bounds for extracting isotopic peaks from the observed spectrum, allowing for some mass inaccuracy specified by mass_accuracy 
extract_isotopes <- function(data, theor_ID_cmpd1, theor_ID_cmpd2, proton_mass, expected_charge_range, deduplicate_fun, mass_accuracy= 50, probability_coverage = 0.99) {
  
  ind1 <- min(which(cumsum(theor_ID_cmpd1$isoDistr)>=probability_coverage))
  ind2 <- min(which(cumsum(theor_ID_cmpd2$isoDistr)>=probability_coverage))
  
  extraction_mass_range <- range(c(theor_ID_cmpd1$masses[1:ind1], theor_ID_cmpd2$masses[1:ind2]))
  
  mz_lb <- (extraction_mass_range[1]-expected_charge_range*proton_mass)/expected_charge_range 
  mz_ub <- (extraction_mass_range[2]-expected_charge_range*proton_mass)/expected_charge_range 
  
  # include mass accuracy in boundaries computations
  mz_lb <- mz_lb - extraction_mass_range[1] * mass_accuracy/10^6
  mz_ub <- mz_ub + extraction_mass_range[2] * mass_accuracy/10^6
  
  charge_df <- data.frame(z = expected_charge_range, mz_lb = mz_lb, mz_ub = mz_ub)
  fun_tmp <- function(x){
    pmap(charge_df, ~{
      #x %>% filter(mz >= ..2, mz <= ..3) %>% mutate(z = ..1, .before = everything())
      x %>% filter(mz >= ..2, mz <= ..3) %>% mutate(z = ..1)
    }) %>% bind_rows() %>% 
      select(z, mz, intensity) # Piotr's update: this order of columns is required
  }
  data <- fun_tmp(data)
  return(list(data = data, charge_df = charge_df))
}
