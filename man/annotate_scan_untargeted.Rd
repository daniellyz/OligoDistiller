% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_scan_untargeted.R
\name{annotate_scan_untargeted}
\alias{annotate_scan_untargeted}
\title{Non-targeted screening from a deconvoluted oligonucleotide spectra}
\usage{
annotate_scan_untargeted(
  scan_processed_aggregated,
  bblock,
  ntheo = 12,
  min_overlap = 0.6,
  max_msigma = 10,
  max_mmw_ppm = 10,
  baseline = 1000
)
}
\arguments{
\item{scan_processed_aggregated}{Data frame representing deconvoluted NMS with the true molecular weight scale. Output of the function process_scan.}

\item{bblock}{Character. Either "DNA" or "RNA". Should reflect the main nucleic acid composition of the strand. Used for monoisotopic peak prediction by Pointless algorithm.}

\item{ntheo}{Integer. Estimated isotope envelop size in number of isotope peaks.}

\item{min_overlap}{Double between 0 and 1. The minimum matching score between experimental and theoretical isotope envelops (known compounds from database/transformation list or unknown predicted by Pointless).}

\item{max_msigma}{Double between 1 and 50. The maximum-allowed deviation between the shapes of experimental and theoretical isotope pattern. Should set higher for noisy or MS/MS data.}

\item{max_mmw_ppm}{Double between 1 and 50. The maximum allowed ppm error between masses in the NMS and theoretical molecular weight of oligonucleotide features. Depend on experimental mass deviation and deconvolution bias.}

\item{baseline}{Numeric. Estimated baseline level (noise) of input spectrum. Depending on instrument and acquisition method. Baseline of MS/MS spectrum is 100 for most instruments.}
}
\description{
The function searches DNA/RNA-like isotope patterns from a deconvoluted oligonucleotide spectra. It provides the monoisotopic molecular weight, average, intensity and envelope likeness of all features detected.
}
\author{
Youzhong Liu, \email{liu-youzhong@hotmail.com}
}
