#-INFO ####
# Function to calculate the REE+Y temperature of coexisting (cogenetic)
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
TREEplcpx <- function(sample = "UNKNOWN",       # sample name | a single character string
                     cpx_m,                     # clinopyroxene major elements | data frame
                     cpx_REEY,                  # clinopyroxene trace elements | data frame
                     pl_m,                      # plagioclase major elements | data frame
                     pl_REEY,                   # plagioclase trace elements | data frame
                     H2O = 0.0,                 # H2O in g/100g
                     P = 2,                     # P in kBar
                     norm_cpx = T,              # when calculating clinopyroxene composition normalize oxides to 100 g/100g | T, F
                     norm_pl = T,               # when calculating plagioclase composition normalize oxides to 100 g/100g | T, F
                     Dcalc = "outer",           # how to handle Dpl-cpx calculation? | "paired" or "outer"
                     REEpresent = 8,            # skip calculation if only less than n REEs+Y are present in the input data | 2 to 15
                     regression = "IWLS",       # type of regression to perform | "IWLS" or "simple"
                     exclude = NULL,            # vector of REEs to exclude from temperature regression | e.g. c("Eu","Y",...) or NULL
                     stripoutlier = F,          # strip REEs from temperature caluclation that classifie as outlier | T, F
                     residual_cutoff = 2,       # values with higher residuals will be filtered
                     inversion_plot = T,        # generate and return the T-inversion plot | T, F
                     REE_plot = T,              # generate normalized REE plot | T, F
                     REE_normalize = "N-MORB"   # Normalize REE to? | N-MORB, C1
                     )  
{  
  #-LOAD PACKAGES & FORMAT BASIC VARIABLES AND SUPPORT FUNCTIONS ####
  require(MASS)
  require(robustbase)
  require(ggplot2)
  require(ggrepel)
  require(EnvStats)
  require(reshape2)
  require(ggpubr)
  
  majors <- c("SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O")
  traces <- c("La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Y","Ho","Er","Tm","Yb","Lu")
  output <- NULL
  
  #----
  #-CHECK INPUT ####
  base::message("** CHECKING DATA ****")
  if (!is.character(sample) & !is.null(sample)){
    stop("***! no sample name given or sample name not a character string.")
  }
  base::message(paste0(
    "**** found ",
    dim(cpx_m)[1],
    " record(s) for clinopyroxene major elements."
  ))
  base::message(paste0(
    "**** found ",
    dim(pl_m)[1],
    " record(s) for plagioclase major elements."
  ))
  base::message(paste0(
    "**** found ",
    dim(cpx_REEY)[1],
    " record(s) for clinopyroxene trace elements."
  ))
  base::message(paste0(
    "**** found ",
    dim(pl_REEY)[1],
    " record(s) for plagioclase trace elements."
  ))
  if (dim(cpx_m)[1] < dim(cpx_REEY)[1] & dim(cpx_m)[1] == 1) {
    base::message(
      paste0(
        "***! clinopyroxene records missmatch! using 1 major record for all trace records."
      )
    )
    cpx_mn <- 1
    cpx_m <- cpx_m[rep(seq_len(nrow(cpx_m)), dim(cpx_REEY)[1]), ]
  } else if (dim(cpx_m)[1] < dim(cpx_REEY)[1] &
             dim(cpx_m)[1] != 1) {
    stop(
      "***! clinopyroxene sample missmatch can't be resolved! \n  ***! clinopyroxene major element samples < clinopyroxene trace element samples \n  ***! either you enter equal number of samples or combine one (major or trace sample) with many "
    )
  }
  
  if (dim(cpx_m)[1] > dim(cpx_REEY)[1] & dim(cpx_REEY)[1] == 1) {
    base::message(
      paste0(
        "***! clinopyroxene records missmatch! using 1 trace record for all major records."
      )
    )
    cpx_REEY <- cpx_REEY[rep(seq_len(nrow(cpx_REEY)), dim(cpx_m)[1]), ]
  } else if (dim(cpx_m)[1] > dim(cpx_REEY)[1] &
             dim(cpx_REEY)[1] != 1) {
    stop(
      "***! clinopyroxene sample missmatch can't be resolved! \n  ***! clinopyroxene major element samples > clinopyroxene trace element samples \n  ***! either you enter equal number of samples or combine one (major or trace sample) with many "
    )
  }

  if (dim(pl_m)[1] < dim(pl_REEY)[1] & dim(pl_m)[1] == 1) {
    base::message(
      paste0(
        "***! plagioclase records missmatch! using 1 major record for all trace records."
      )
    )
    pl_mn <- 1
    pl_m <- pl_m[rep(seq_len(nrow(pl_m)), dim(pl_REEY)[1]), ]
  } else if (dim(pl_m)[1] < dim(pl_REEY)[1] & dim(pl_m)[1] != 1) {
    stop(
      "***! plagioclase sample missmatch can't be resolved! \n  ***! plagioclase major element samples < plagioclase trace element samples \n  ***! either you enter equal number of samples or combine one (major or trace sample) with many "
    )
  }
  if (dim(pl_m)[1] > dim(pl_REEY)[1] & dim(pl_REEY)[1] == 1) {
    base::message(
      paste0(
        "***! plagioclase records missmatch! using 1 trace record for all major records."
      )
    )
    pl_REEY <- pl_REEY[rep(seq_len(nrow(pl_REEY)), dim(pl_m)[1]), ]
  } else if (dim(pl_m)[1] > dim(pl_REEY)[1] &
             dim(pl_REEY)[1] != 1) {
    stop(
      "***! plagioclase sample missmatch can't be resolved! \n  ***! plagioclase major element samples > plagioclase trace element samples \n  ***! either you enter equal number of samples or combine one (major or trace sample) with many "
    )
  }
  if (dim(pl_m)[1] != dim(cpx_m)[1] & Dcalc != "outer") {
    Dcalc <- "outer"
    message(
      "***! paired calculation not possible with diferent number of records for clinopyroxene and plagioclase, switching to 'outer'."
    )
  }
  message("** SET OPTIONS ****")
  if (base::is.null(sample)) {
    sample <- "unknown"
    base::message(paste0("***! you have not entered a sample name, setting it to 'unknown'."))
  } else {
    base::message(paste0("**** sample = ", sample))
  }
  message(paste0("**** H2O = ", H2O, " g/100g"))
  message(paste0("**** p = ", P, " kBar"))
  message(paste0("**** norm_cpx = ", norm_cpx))
  message(paste0("**** norm_pl = ", norm_cpx))
  base::message(paste0("**** Dcalc = ", Dcalc))
  if (Dcalc == "paired") {
    CalcMode <- rep(Dcalc == "paired", dim(cpx_m)[1])
  } else if (Dcalc == "outer") {
    CalcMode <- rep(!Dcalc == "outer", dim(cpx_m)[1] * dim(pl_m)[1])
  } else {
    stop("***! Couldn't find mode for Dcalc, set it to 'paired' or 'outer'")
  }
  if (length(exclude) != 0){
    message(paste0("***! adjusting the minimum number of REEs present for calculation to account for the excluded REEs"))
    REEpresent <- REEpresent-length(exclude)
  }
  message(paste0("**** REEpresent = ", REEpresent))
  message(paste0("**** exclude = ", paste(exclude, collapse=" ")))
  message(paste0("**** stripoutlier = ", stripoutlier))
  message(paste0("**** inversion_plot = ", inversion_plot))
  
  message("** STARTING CALCULATION ****")
  #----
  #-ATOM WEIGHTS AND STOICHIOMETRY ####
  # see: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl "Standard Atomic Weights"
  mol_oxygen <- 15.999
  mol_cat_an <- data.frame(
    SiO2 = c(28.084 + mol_oxygen * 2, 1, 2),
    TiO2 = c(47.867 + mol_oxygen * 2, 1, 2),
    Al2O3 = c(26.981 * 2 + mol_oxygen * 3, 2, 3),
    Cr2O3 = c(51.996 * 2 + mol_oxygen * 3, 2, 3),
    FeO = c(55.845 + mol_oxygen, 1, 1),
    MnO = c(54.938 + mol_oxygen, 1, 1),
    MgO = c(24.304 + mol_oxygen, 1, 1),
    CaO = c(40.078 + mol_oxygen, 1, 1),
    Na2O = c(22.989 * 2 + mol_oxygen, 2, 1),
    K2O = c(39.098 * 2 + mol_oxygen, 2, 1)
  )
  mol_cat_an_cpx <- mol_cat_an[majors %in% names(cpx_m)]
  mol_cat_an_pl <- mol_cat_an[majors %in% names(pl_m)]
  #----
  #-RECALULATE MINERAL FORMULA FOR CPX ####
  message("**** recalculating cpx mineral formula.")
  cpx <- cpx_m[names(cpx_m) %in% majors]
  cpx[is.na(cpx)] <- 0
  cpx_t <- rowSums(cpx)
  if (norm_cpx == T) {cpx <- cpx * 100 / cpx_t}
  o <- cpx / mol_cat_an_cpx[1, ][rep(seq_len(nrow(mol_cat_an_cpx[1, ])), dim(cpx_m)[1]), ] * mol_cat_an_cpx[3, ][rep(seq_len(nrow(mol_cat_an_cpx[3, ])), dim(cpx_m)[1]), ]
  o_f <- 6 / rowSums(o)
  cat_cpx <- o * o_f / mol_cat_an_cpx[3, ][rep(seq_len(nrow(mol_cat_an_cpx[3, ])), dim(cpx_m)[1]), ] * mol_cat_an_cpx[2, ][rep(seq_len(nrow(mol_cat_an_cpx[2, ])), dim(cpx_m)[1]), ]
  cpx[["total"]] <- rowSums(cpx)
  cpx[["AlIV"]] <-  ifelse((cat_cpx$SiO2 < 2) & (cat_cpx$Al2O3 > 2 - cat_cpx$SiO2),2 - cat_cpx$SiO2,ifelse(cat_cpx$Al2O3 <= 2 - cat_cpx$SiO2,cat_cpx$Al2O3,stop("***! can't resolve cpx formula, check if you hava a stoichiometric cpx, check if you assigned your variables the right way, cahnge norm_cpx.")))
  cpx[["AlVI"]] <-  ifelse((cat_cpx$SiO2 < 2) & (cat_cpx$Al2O3 > 2 - cat_cpx$SiO2),cat_cpx$Al2O3 - cpx$AlIV,ifelse(cat_cpx$Al2O3 <= 2 - cat_cpx$SiO2,0,stop("***! can't resolve cpx formula, check if you hava a stoichiometric cpx, check if you assigned your variables the right way, cahnge norm_cpx.")))
  cpx["MgM2"] <- (1 - (cat_cpx$CaO + cat_cpx$Na2O + cat_cpx$K2O + cat_cpx$MnO)) / (1 + cat_cpx$FeO / cat_cpx$MgO)
  cpx[["MgNo"]] <-  cat_cpx$MgO / (cat_cpx$FeO + cat_cpx$MgO) * 100
  #----
  #-RECALULATE MINERAL FORMULA FOR PL ####
  message("**** recalculating pl mineral formula.")
  pl <- pl_m[names(pl_m) %in% majors]
  pl[is.na(pl)] <- 0
  pl_t <- rowSums(pl)
  if (norm_pl == T) {pl <- pl * 100 / pl_t}
  o <- pl / mol_cat_an_pl[1,][rep(seq_len(nrow(mol_cat_an_pl[1,])), dim(pl_m)[1]), ] * mol_cat_an_pl[3,][rep(seq_len(nrow(mol_cat_an_pl[3,])), dim(pl_m)[1]), ]  
  o_f <- 8 / rowSums(o)
  cat_pl <- o * o_f / mol_cat_an_pl[3,][rep(seq_len(nrow(mol_cat_an_pl[3,])), dim(pl_m)[1]), ] * mol_cat_an_pl[2,][rep(seq_len(nrow(mol_cat_an_pl[2,])), dim(pl_m)[1]), ]
  pl[["total"]] <- rowSums(pl)
  pl[["XCa"]] <- cat_pl$CaO
  pl[["XNa"]] <- cat_pl$Na2O
  pl[["XK"]] <- cat_pl$K2O
  pl[["An"]] <- cat_pl$CaO / (cat_pl$CaO + cat_pl$Na2O + cat_pl$K2O) * 100
  #----
  #-CALCULATE TEMPERATURE ####
  message("**** calculating inversion indices.")
  # see: http://abulafia.mt.ic.ac.uk/shannon/ptable.php "Charge = 3, Coordination = VIII" 
  irj <-
    c(
      La = 1.16,
      Ce = 1.143,
      Pr = 1.126,
      Nd = 1.109,
      Sm = 1.079,
      Eu = 1.066,
      Gd = 1.053,
      Tb = 1.04,
      Dy = 1.027,
      Y = 1.019,
      Ho = 1.015,
      Er = 1.004,
      Tm = 0.994,
      Yb = 0.985,
      Lu = 0.977
    )
  
  pl_REEY <- pl_REEY[names(pl_REEY) %in% traces]
  cpx_REEY <- cpx_REEY[names(cpx_REEY) %in% traces]
  #----
  #--X_melt_H2O CALCULATION ####
  # following: Wood, B. J., & Blundy, J. D. (2002), The effect of H2O on crystal-melt partitioning of trace elements, Geochimica Et Cosmochimica Acta, 66(20), 3647–3656, http://doi.org/10.1016/S0016-7037(02)00935-3
  X_melt_H2O <- ifelse(H2O > 0,
                       1.70685747186617E-04 * H2O ^ 3 - 6.831778563526E-03 * H2O ^ 2 + 1.09999984719062E-01 * H2O + 9.71902247238525E-04,
                       0)
  #----
  #--ld(Dj)-A CALCULATION ####
  A <- 1.60468782786597E+01 + 7.13615731307051
  Ao <- -5.17178970845966 * pl$XCa * pl$XCa
  Ac <- 4.3709642058855 * cpx$AlIV + 1.98130564180907 * cpx$MgM2 - 9.08067851952162E-01 * X_melt_H2O
  A <- ifelse(CalcMode,
              A + Ao - Ac,
              A + as.vector(outer(as.numeric(Ao), as.numeric(Ac), function(X, Y) X - Y))
              )
  xREE <- NULL
  for (i in traces){
    xREE[[i]] <- ifelse(CalcMode,
                        (log(pl_REEY[[i]]/cpx_REEY[[i]]) - A),
                        (log(as.vector(outer(as.numeric(pl_REEY[[i]]),as.numeric(cpx_REEY[[i]]), function(X, Y) X / Y))) - A)
                        )
  }
  xREE <- as.data.frame(xREE)
  xREEin <- as.data.frame(xREE)
  xREEout <- as.data.frame(xREE)
  xREEin[exclude] <- NA
  xREEex <- xREE[names(xREE) %in% exclude]
  #----
  #--Bj/1000 CALCULATION ####
  PI <- 3.14159265358979
  Avogadro <- 602.214076
  r0_pl <- rep(1.17908321870255,dim(pl_m)[1]) #Å ionic radius of a strain-free cation in pl
  E_pl <- 1.96175540567148E+02 #GPa Young's modulud for pl 
  r0_cpx <- 1.06596147404672 - 1.03654121861864E-01 * cpx$AlVI - 2.1180375646043E-01 * cpx$MgM2
  E_cpx <- -1.99606952132459E+03 + 2.2719097731105E+03 * r0_cpx
  if (P > 0){
    PGPa <- P / 10
  }
  yREE <- NULL
  B_cpx <- NULL
  B_pl <- NULL
  for (i in traces){
    B_pl[[i]] <- 4 * PI * E_pl * Avogadro * (r0_pl / 2 * (r0_pl - irj[[i]]) ^ 2 - (r0_pl - irj[[i]]) ^ 3 / 3)
    B_cpx[[i]] <- 4 * PI * E_cpx * Avogadro * (r0_cpx / 2 * (r0_cpx - irj[[i]]) ^ 2 - (r0_cpx - irj[[i]]) ^ 3 / 3)
    yREE[[i]] <- 
      ifelse(CalcMode,
             (((-1.94543445294503E+05 - 7.18648736931024E+04 - 1.17036448671446E+04 * PGPa^2  - B_pl[[i]]+ B_cpx[[i]]) / 8.3145) / 1000),
             (((-1.94543445294503E+05 - 7.18648736931024E+04 - 1.17036448671446E+04 * PGPa^2  - as.vector(outer(as.numeric(+B_pl[[i]]), as.numeric(-B_cpx[[i]]), function(X, Y) X + Y))) / 8.3145) / 1000)
             )
  }

  yREE <- as.data.frame(yREE)
  yREEin <- as.data.frame(yREE)
  yREEout <- as.data.frame(yREE)
  yREEin[exclude] <- NA
  yREEex <- yREE[names(yREE) %in% exclude]
  #----
  #--STRIP OUTLIERS ####
  TREEplcpx_T <- NULL
  TREEplcpx_err <- NULL
  TREEplcpx_rot <- NULL
  if (stripoutlier == T){
      message("**** cleaning data for outliers by applying a Rosner's test to the rotation residuals")
    analytes_in <- names(xREEin)
    n_col <- dim(xREEin)[2]
    n_row <- dim(xREEin)[1]
    yREEout <- unlist(yREEout)
    yREEin <- unlist(yREEin)
    xREEout <- unlist(xREEout)
    xREEin<- unlist(xREEin)
    
    resi <-  residuals(rlm(yREEin~xREEin, method = "MM", scale.est = "MAD",psi = psi.huber,na.action="na.exclude"))^2

    xREEin[resi > residual_cutoff] <- NA
    yREEin[resi > residual_cutoff] <- NA
    
    xREEout[resi < residual_cutoff] <- NA
    yREEout[resi < residual_cutoff] <- NA
    
    xREEin <- as.data.frame(matrix(xREEin, nrow = n_row, ncol = n_col))
    yREEin <- as.data.frame(matrix(yREEin, nrow = n_row, ncol = n_col))
    names(xREEin) <- analytes_in
    names(yREEin) <- analytes_in
    
    xREEout <- as.data.frame(matrix(xREEout, nrow = n_row, ncol = n_col))
    yREEout <- as.data.frame(matrix(yREEout, nrow = n_row, ncol = n_col))
    names(xREEout) <- analytes_in
    names(yREEout) <- analytes_in

    message(paste0(capture.output(colSums(is.na(xREEin))), collapse = "\n"))
  }
  #--INVERSION REGRESSION #### 
  message("**** starting temperature regression.")
  n = 0
  for (i in 1:dim(xREEin)[1]) {
    if (rowSums(!is.na(xREEin[i,]))[[1]] < REEpresent | rowSums(!is.na(yREEin[i,]))[[1]] < REEpresent){
      next()
    }
    n = n+1
    x <- as.numeric(xREEin[i,])
    y <- as.numeric(yREEin[i,])
    
    if(regression == "IWLS"){
      fixed_model <- MASS::rlm(unlist(y)~unlist(x)+0, method = "MM", scale.est = "MAD", psi = psi.huber)
      loose_model <- MASS::rlm(unlist(y)~unlist(x), method = "MM", scale.est = "MAD",psi = psi.huber)
      
      } else {
        fixed_model <- stats::lm(unlist(y)~unlist(x)+0)
        loose_model <- stats::lm(unlist(y)~unlist(x))
        }
    TREEplcpx_T[[i]] <-  summary(fixed_model)$coefficients[1]
    TREEplcpx_err[[i]] <-  summary(fixed_model)$coefficients[2]
    TREEplcpx_rot[[i]] <- summary(loose_model)$coefficients[1]
  }
  output[["Temperature"]] <- unlist(TREEplcpx_T)*1000-273.15
  output[["median_Temperature"]] <- stats::median(unlist(TREEplcpx_T), na.rm = T)*1000-273.15
  output[["median_model_error"]] <- stats::median(unlist(TREEplcpx_err), na.rm = T)*1000
  output[["sample_Temperature_range"]] <- (base::max(unlist(TREEplcpx_T), na.rm = T)-base::min(unlist(TREEplcpx_T), na.rm = T))/2*1000
  output[["median_rotation"]] <- stats::median(unlist(TREEplcpx_rot), na.rm = T)
  output[["Bj1000"]] <- yREEin
  output[["lnDjA"]] <- xREEin
  output[["MgNo"]] <- cpx$MgNo
  output[["An"]] <- pl$An
  output[["xCa"]] <- pl$XCa
  output[["cpxAlIV"]] <- cpx$AlIV
  #----
  #-CREATE PLOTS ####
  message("**** generating inversion plot.")
  df <- data.frame(element = reshape2::melt(xREEin)$variable, x_in = reshape2::melt(xREEin)$value, y_in = reshape2::melt(yREEin)$value)
  if (stripoutlier == T){
    df <- base::cbind(df,x_out = reshape2::melt(xREEout)$value, y_out = reshape2::melt(yREEout)$value)
  }
  if (inversion_plot == T){
      if (stripoutlier == T){
        p_inv <- 
          ggplot()+
          geom_point(data = df, aes(x = x_in, y = y_in), size = 2, color = "black", shape = 1, stroke = 0.2, alpha = 1)+
          geom_point(data = df, aes(x = x_out, y = y_out), size = 2, color = "red", shape = 4, stroke = 0.2, alpha = 1)+
          geom_point(aes(x = robustbase::colMedians(as.matrix(xREEin), na.rm = T), y = robustbase::colMedians(as.matrix(yREEin),na.rm = T)), size = 3, color = "black", fill = "yellow", shape = 23, alpha = 1)
      } else {
        p_inv <- 
          ggplot()+
          geom_point(data = df, aes(x = x_in, y = y_in), size = 2, color = "black", shape = 1, stroke = 0.2, alpha = 1)+
          geom_point(aes(x = robustbase::colMedians(as.matrix(xREEin), na.rm = T), y = robustbase::colMedians(as.matrix(yREEin),na.rm = T)), size = 3, color = "black", fill = "yellow", shape = 23, alpha = 1)
      }
    
      if(regression == "IWLS"){
      p_inv$layers <- c(
        geom_smooth(aes(x = reshape2::melt(robustbase::colMedians(as.matrix(xREEin),na.rm = T))$value, y = reshape2::melt(robustbase::colMedians(as.matrix(yREEin),na.rm = T))$value),
                    method = MASS::rlm,
                    formula = y ~ x + 0,
                    method.args=c(scale.est = "MAD",method = "MM", psi = psi.huber),
                    fullrange = T, se = T, fill="lightblue", color = "black", size = 1, level = 0.95, na.rm = T),
        geom_smooth(aes(x = reshape2::melt(robustbase::colMedians(as.matrix(xREEin),na.rm = T))$value, y = reshape2::melt(robustbase::colMedians(as.matrix(yREEin),na.rm = T))$value),
                    method = MASS::rlm,
                    formula = y ~ x,
                    method.args=c(scale.est = "MAD",method = "MM", psi = psi.huber),
                    fullrange = T, se = F, linetype = "dashed", color = "black", size = 0.4, na.rm = T),
        p_inv$layers)
      } else {
        p_inv$layers <- c(
        geom_smooth(data = df, aes(x = x_in,y = y_in), method = stats::lm, formula = y ~ x + 0, fullrange = T, se = T, fill="lightblue", color = "black", size = 1, level = 0.95, na.rm = T),
        geom_smooth(data = df, aes(x = x_in,y = y_in), method = stats::lm, formula = y ~ x, fullrange = T, se = F, color = "black",linetype = "dashed",size = 1, level = 0.95, na.rm = T),
        p_inv$layers)
      }
    if(!is.null(exclude)){
      p_inv$layers <- c(
        p_inv$layers,
        geom_point(aes(x = reshape2::melt(as.matrix(xREEex))$value, y = reshape2::melt(as.matrix(yREEex))$value), size = 2, color = "red", shape = 4, stroke = 0.2, alpha = 1)
      )
    }
     p_inv <-  p_inv+
    theme_classic()+
      theme(axis.title.x = element_text(color = "darkblue", size = 10, face = "bold"),
            axis.title.y = element_text(color = "darkblue", size = 10, face = "bold"),
            axis.text = element_text(color = "black", size = 9, face = "bold"),
            panel.grid.major = element_blank(),
            plot.title = element_text(
              color = "darkblue",
              size = 14,
              face = "bold",
              hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1)
      )+
      ylab(paste0("Bj / 1000"))+
      xlab(paste0("ln(Dj) - A"))+
      scale_x_continuous(expand=c(0,0), limits=c(-60,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(-60,0)) +
      coord_cartesian(ylim = c(-36, -28), xlim=c(-28,-17))+
      ggplot2::annotate("label_repel", x = -27.5, y = -28, fill = "white", label = paste0(
                                                                               "sample: ", sample,
                                                                               "\n md model: ", signif(output$median_Temperature,4),
                                                                               " °C\n se model: ± ",signif(output$median_model_error,2),
                                                                               " °C\n range sample: ± ",signif(output$sample_Temperature_range,2),
                                                                               " °C\n ",Dcalc, " n: ",n,
                                                                               "\n -------------------- \n Thy09: ",signif(899+3.6*stats::median(output$An, na.rm = T),4),
                                                                               " °C\n Morse10: ",signif(865+4.6*stats::median(output$An, na.rm = T),4),
                                                                               " °C\n V&M13: ",signif(6.8*(stats::median(output$An, na.rm = T)-75)+1150,4),
                                                                               " °C\n REE (An): ",signif(706+6.1*stats::median(output$An, na.rm = T),4),
                                                                               " °C\n -------------------- \n pl An: ",signif(stats::median(output$An,na.rm = T),3),
                                                                               " \n cpx Mg#: ", signif(stats::median(output$MgNo, na.rm = T),3),
                                                                               " \n H20: ",H2O," [wt%]\n p: ",signif(P,2),
                                                                               " [kBar]"),
               hjust = 0.5,vjust= 0.5, segment.size = 0, size = 3)+
      geom_label_repel(data = NULL, inherit.aes = F, aes(x = robustbase::colMedians(as.matrix(xREE), na.rm = T)[c(1,5,6,10,length(yREE))], y = robustbase::colMedians(as.matrix(yREE), na.rm = T)[c(1,5,6,10,length(yREE))]), label = colnames(yREE)[c(1,5,6,10,length(yREE))], label.padding = 0.2, box.padding = 0.2, segment.size = 0.4, force = 2, size = 2.2 )
    } 
  if (REE_plot == T){
    if (REE_normalize == "N-MORB"){
    normal <- c(La = 4.19, Ce = 12.42, Pr = 1.98, Nd = 10.66, Sm = 3.48,Eu = 1.26,Gd = 4.55,Tb = 0.82,Dy = 5.5, Y = 33.2,Ho = 1.18,Er = 3.42,Tm = 0.52,Yb = 3.28,Lu = 0.48)
    } else if (REE_normalize == "C1"){
    normal <- c(La = 0.237, Ce = 0.613, Pr = 0.0928, Nd = 0.457, Sm = 0.148, Eu = 0.0563, Gd = 0.199, Tb = 0.0361, Dy = 0.246, Y = 1.57, Ho = 0.0546, Er = 0.16, Tm = 0.0247, Yb = 0.161,Lu = 0.0246)
    }
    x_order <-   c("La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Y","Ho","Er","Tm","Yb","Lu")
      df_cpx           <- reshape2::melt(sweep(t(as.data.frame(t(cpx_REEY),row.names = x_order)), 2, normal, "/"))
      names(df_cpx)    <- c("Sample","Element","Value")
      df_cpx           <- as.data.frame(df_cpx)

      df_pl           <- reshape2::melt(sweep(t(as.data.frame(t(pl_REEY),row.names = x_order)), 2, normal, "/"))
      names(df_pl)    <- c("Sample","Element","Value")
      df_pl         <- as.data.frame(df_pl)
      
      p_REE <- 
      ggplot()+
        geom_line(
          data = df_cpx,
          aes(
            x = factor(Element, levels = x_order),
            y = Value,
            group = Sample
          ),
          linetype = "solid",
          size = 0.3,
          color = "gray"
        )+
        geom_line(
          data = df_pl,
          aes(
            x = factor(Element, levels = x_order),
            y = Value,
            group = Sample
          ),
          linetype = "solid",
          size = 0.3,
          color = "gray"
        )+
        geom_line(aes(x = factor(x_order, levels = names(colMedians(as.matrix(cpx_REEY),na.rm = T))),
                      y = colMedians(as.matrix(cpx_REEY), na.rm = T)/normal,
                      group = 1))+
        geom_line(aes(x = factor(x_order,levels = names(colMedians(as.matrix(pl_REEY)))),
                      y = colMedians(as.matrix(pl_REEY),na.rm = T)/normal,
                      group = 1))+
        geom_point(aes(x = factor(x_order, levels = names(colMedians(as.matrix(cpx_REEY),na.rm = T))),
                       y = colMedians(as.matrix(cpx_REEY), na.rm = T)/normal, shape = "cpx"),color = "black", fill = "yellow", size = 3)+
        geom_point(aes(x = factor(x_order, levels = names(colMedians(as.matrix(pl_REEY)))),
                       y = colMedians(as.matrix(pl_REEY),na.rm = T)/normal, shape = "pl"),color = "black", fill = "yellow", size = 3)+
        scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10 ^ .x)))+
        #scale_fill_manual(labels = c("pl", "cpx"), values = c("pl" = "yellow","cpx" = "yellow")) +
        #scale_color_manual(values = c("pl" = "black", "cpx" = "black")) +
        scale_shape_manual(name = NULL, breaks = c("pl","cpx"), values = c("pl" = 24, "cpx" = 21))+
        annotation_logticks(sides = "l")+
        theme_classic() +
        theme(
          plot.title = element_text(
            color = "black",
            size = 14,
            face = "bold",
            hjust = 0.5
          ),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "darkblue", size = 10, face = "bold"),
          axis.text = element_text(color = "black", size = 9, face = "bold"),
          panel.grid.major.y = element_line(colour = "gray22",
                                            size = 0.1),
          panel.border = element_rect(
            colour = "black",
            fill = NA,
            size = 1
          ),
          legend.position = c(0.08, 0.1),
          legend.background = element_rect(
                                           size=0.5, linetype="solid", 
                                           colour ="black")
        ) +
        ylab(paste0(REE_normalize, " normalized REE mass fraction [µg/g]"))
  }
  if (REE_plot == T & inversion_plot == T){
    p <- ggarrange(p_REE, p_inv, widths = c(3, 2), heights = c(1,1),
              ncol = 2, nrow = 1, align = "h")
    output[["inversion_plot"]] <- p
  } else if (REE_plot == F & inversion_plot == T){
    output[["inversion_plot"]] <- p_inv
  } else {
    message("***! no inversion plot in output")
  }
 
  
  
  #----
  message("** DONE ****")
  suppressWarnings(return(output))
}
