#-INFO ####
# Example of using the function to calculate the REE+Y temperature of coexisting (cogenetic)
# cpx-pl pairs or an outer matrix of all data collected for one sample
# 
# based on:
# Sun, C., & Liang, Y. (2017).  
# A REE-in-plagioclase–clinopyroxene thermometer for crustal rocks
# CMP, 172(4), 1234.
# http://doi.org/10.1007/s00410-016-1326-9
#
# written by:
# Samuel Müller
# Institute of Geosceiences, Kiel University
# samuel.mueller@ifg.uni-kiel.de
# 
# published in (cite as): 
# Müller, S., Garbe-Schönberg, D., Koepke, J., Hoernle, K. (2021).
# A reference section through fast-spread lower oceanic crust, Wadi Gideah, Samail Oophiolite (Sultanate of Oman):
# Trace Element Systematics and Crystallization Temperatures – implications for hybrid crustal accretion
# Journal of Geophysical Research: Solid Earth xx(xx), xx.
# DOI
#
# Version 0.3
# for help and bugs don't hazitate to contact me via mail
#----
# get the current path of the example.R file ! this only works in R-Studio
rm(list=ls())
userpath <- dirname(rstudioapi::getSourceEditorContext()$path)

# source the function
source(paste0(userpath,"/TREEplcpx.R"))

# read the sample data
cpx_major <- read.csv(paste0(userpath,"/cpx_m.csv"))
cpx_trace <- read.csv(paste0(userpath,"/cpx_REEY.csv"))
pl_major <- read.csv(paste0(userpath,"/pl_m.csv"))
pl_trace <- read.csv(paste0(userpath,"/pl_REEY.csv"))

# calculate temperatures

TREEplcpx(
  sample = "OM12-025",
  cpx_m = cpx_major,
  cpx_REEY = cpx_trace,
  pl_m = pl_major,
  pl_REEY = pl_trace,
  H2O = 1,
  P = 2,
  norm_cpx = T,
  norm_pl = T,
  Dcalc = "outer",
  REEpresent = 8,
  regression = "IWLS",
  exclude = NULL,
  stripoutlier = T,
  residual_cutoff = 2,
  inversion_plot = T,
  REE_plot = T,
  REE_normalize = "N-MORB"
)
