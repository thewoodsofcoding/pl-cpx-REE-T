# pl-cpx-REE-T
R function to calculate plagioclase-clinopyroxene REE exchange temperatures in a flexible way.
For details see: Müller, S., Garbe-Schönberg, D., Koepke, J., Hoernle, K. (202x). A reference section through fast-spread lower oceanic crust, Wadi Gideah, Samail ophiolite (Sultanate of Oman): Trace Element Systematics and Crystallization Temperatures – implications for hybrid crustal accretion. Journal of Geophysical Research-Solid Earth, xxx(x), DOI

![](https://raw.githubusercontent.com/thewoodsofcoding/pl-cpx-REE-T/main/img/inversion.jpg)

The calculations are base on: Sun, C., & Liang, Y. (2017). A REE-in-plagioclase–clinopyroxene thermometer for crustal rocks. Contributions to Mineralogy and Petrology, 172(4), 24. [DOI](http://doi.org/10.1007/s00410-016-1326-9)

# getting started
1. Download and isntall [R-Studio](https://www.rstudio.com/products/rstudio/download/)
2. Download this repository as zip file and extract it.
3. Open example.R file in R-Studio and hit Source.
4. You have to make sure that all required packages are installed, if this is not the case type the following into the R-Studio console:
```
install.packages("MASS","robustbase","ggplot2","ggrepel","EnvStats","reshape2","ggpubr")
```
5. Beyond the example.R, you can call the R function to calculate REE temperatures with `TREEplcpx()` after you have sourced it to your working environment via `source(paste0(USERPATH,"/TREEplcpx.R"))`. Replace the USERPATH with the path to the TREEplcpx.R file of the downloaded zip folder.

# options
| option          | values           |default    | description                                                                  
|:----------------|:-----------------|:----------|:-----------------------------------------------------------------------------|
| sample          | character string | UNKNOWN   | sample name                                                                  
| cpx_m           | data frame       |           | clinopyroxene major elements                                                 
| cpx_REEY        | data frame       |           | clinopyroxene trace elements                                                 
| pl_m            | data frame       |           | plagioclase major elements                                                   
| pl_REEY         | data frame       |           | plagioclase trace elements                                                   
| H2O             | 0 to 100         | 0.0       | H2O in g/100g                                                                
| P               | > 0              | 2         | P in kBar                                                                    
| norm_cpx        | TRUE, FALSE      | TRUE      | when calculating clinopyroxene composition normalize oxides to 100 g/100g?   
| norm_pl         | TRUE, FALSE      | TRUE      | when calculating plagioclase composition normalize oxides to 100 g/100g?     
| Dcalc           | paired, outer    | outer     | Dpl-cpx calculation? Calcualte pairwise or generate an outer matrix?         
| REEpresent      | 2 to 15          | 8         | skip calculation if less than n REEs+Y are present in the input data.         
| regression      | simple, IWLS     | IWLS      | type of regression to perform regular (simple) or robust (IWLS) regression?
| exclude         | REE              | NULL      | vector of REEs to exclude from temperature regression e.g. c("Eu","Y").
| stripoutlier    | TRUE, FALSE      | FALSE     | strip data from temperature caluclation that classifie as outlier?
| residual_cutoff | > 1              | 2         | values with higher residuals will be outliers.
| inversion_plot  | TRUE, FALSE      | TRUE      | generate and return the T-inversion plot?
| REE_plot        | TRUE, FALSE      | TRUE      | generate and return normalized REE plot?
| REE_normalize   | N-MORB, C1       | N-MORB    | normalize REE to N-MORB ([Gale et al., 2013](http://doi.org/10.1029/2012GC004334)) or C1 ([McDonough & Sun, 1995](http://doi.org/10.1016/0009-2541(94)00140-4))
