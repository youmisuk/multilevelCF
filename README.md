# Random Forests Approach for Causal Inference with Clustered Observational Data

Youmi Suk, Hyunseung Kang, and Jee-Seon Kim

## Overview

This paper investigates using one particular machine learning (ML) method based on random forests known as Causal Forests to estimate
treatment effects in multilevel observational data. We conduct simulation studies under different types of multilevel data, including two-level, three-level, and cross-classified data. Our simulation study shows that when the ML method is supplemented with estimated
propensity scores from multilevel models that account for clustered/hierarchical structure, the modified ML method outperforms pre-existing methods in a wide variety of settings. We demonstrate our findings by examining the effect of private math lessons in the Trends in International Mathematics and Science Study (TIMSS) data, a large-scale educational assessment where students are nested within schools.

Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis using TIMSS data. 

## Simulation Study

* `DGM_ClusteredData.R`  

   This `R` file includes data generating codes for two-level data, three-level data, and cross-classified data.
 
```R
twolevel.pop()
threelevel.pop()
ccrem.pop() 
```

* `Simu_MultilevelCF.R`
 
   This `R` file includes simulation codes with different multilevel data structures. One can change simulation parameters: `smpl.size`, `sd.c`, `tau1`, and `tau2`.  For more details on simulation condtions, see our paper, https://psyarxiv.com/xgq2k/.


## TIMSS Data Study

* `TIMSS2015Korea_Math_complete.csv`

  This is our complete data. The original 2015 TIMSS dataset is available at https://timssandpirls.bc.edu/timss2015/international-database/. 

* `TIMSS2015Korea_Math_DataAnalysis.R` 
 
   This `R` file can be used to replicate our data analysis.
