#' Prediction ESI fragments
#'
#' The function creates fragment ions for a sequence given
#' 
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#' 
#' @importFrom OrgMassSpecR MonoisotopicMass
#' @export

predict_esi_frag<-function(test_seq = "OH-Ad-Ad-Ad-Ad-Ad-Ad-OH"){
  
  # Symbols of subunits
  
  ref_unit2 = read.csv("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/RNA-DNA2.txt", sep = "\t")
  sugar2 = lapply(ref_unit2$sugar, function(x) unlist(ListFormula1(x)))
  phosphate2 = lapply(ref_unit2$phosphate, function(x) unlist(ListFormula1(x)))
  base2 = lapply(ref_unit2$base, function(x) unlist(ListFormula1(x)))
  
  water = unlist(ListFormula1("H2O"))
  oh = unlist(ListFormula1("HO"))
  proton = unlist(ListFormula1("H1"))
  PO2 = unlist(ListFormula1("PO2"))
  
  # From sequence to FLP formula:
  
  units = strsplit(test_seq,"-")[[1]]
  ustart = units[1]
  ustart = unlist(ListFormula1(ustart))
  uend = units[length(units)]
  uend = unlist(ListFormula1(uend))
  
  test_seq = units[2:(length(units)-1)]  # Sequence without heads
  NT = length(test_seq)
  
  vvv = match(test_seq, ref_unit2$alias)
  if (NA %in% vvv){stop("Unknown unit detected!")}
  tmp_ref = ref_unit2[vvv,,drop=FALSE]

  tmp_P2 = rep(0, length(ustart))
  tmp_P2 =  Reduce("+", phosphate2[vvv[1:(NT-1)]])
  flp = Reduce("+", sugar2[vvv]) + Reduce("+", base2[vvv]) + ustart + uend + tmp_P2
  flp = flp - water*(NT-1) - 2*oh 
  flp_formula = generate_formula(flp)

  # Generate formula:
  
  fragment_ions_labels = c("FLP")
  fragment_formula = c(flp_formula)
  
  for (i in 1:(NT-1)){
    
    vvv = match(test_seq[1:i], ref_unit2$alias)
    tmp_ref = ref_unit2[vvv,,drop=FALSE]
    
    al = paste0("a", i)
    tmp_P2 = rep(0, length(ustart))
    if (i>1){tmp_P2 =  Reduce("+", phosphate2[vvv[1:(i-1)]])}
    alf = Reduce("+", sugar2[vvv]) + Reduce("+", base2[vvv]) + ustart + tmp_P2
    alf = alf - water*i - oh
    af = generate_formula(alf)

    aBl = paste0("aB", i)
    aBlf = alf - base2[[vvv[i]]]
    aBf = generate_formula(aBlf)
    
    bl = paste0("b", i)
    blf = alf + water
    bf = generate_formula(blf)
    
    cl = paste0("c", i)
    clf = blf + PO2 - proton
    cf = generate_formula(clf)
    
    dl = paste0("d", i)
    dlf = clf + water
    df = generate_formula(dlf)
    
    wl = paste0("w", NT - i)
    wf = generate_formula(flp - alf)
    
    xl = paste0("x", NT - i)
    xf = generate_formula(flp - blf)
    
    yl = paste0("y", NT - i)
    yf = generate_formula(flp - clf)
    
    zl = paste0("z", NT - i)
    zf = generate_formula(flp - dlf)
    
    fragment_ions_labels = c(fragment_ions_labels, al, bl, cl, dl, 
                             aBl, wl, xl, yl, zl)
    
    fragment_formula = c(fragment_formula, af, bf, cf, df, 
                         aBf, wf, xf, yf, zf)
  }
  
  # Create formula table:
  
  summary_formula = cbind.data.frame(ID = 1:length(fragment_ions_labels),
                                     CPD = fragment_ions_labels, FORMULA = fragment_formula)
  
  summary_formula$NM = sapply(fragment_formula, function(x) MonoisotopicMass(ListFormula1(x)))
  summary_formula$NEG = summary_formula$NM - 1.00726 
  
  return(summary_formula)
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

generate_formula<-function(x){
  elements = names(x)
  valid = which(x>0)
  elements = elements[valid]
  x = as.numeric(x[valid])
  y = paste0(paste0(elements,x),collapse="")
  return(y)
}
