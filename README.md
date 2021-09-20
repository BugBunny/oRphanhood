---
title: "oRphanhood"
author: "Ian Timæus"
date: "20/09/2021"
output: html_document
---

The R script calculates the proportions of people with living mothers or fathers by age in a set of stable populations defined by 2-parameter relational logit model life tables, 2-parameter relational Gompertz fertility models, and their growth rates, *r*. These populations are then used to model and predict:

1.  the relationship between life table survivorship (*lm*+*n*/*lm*, where *m* is 25 for women and 35 for men) and the proportion of the population by five-year age group with a living father or mother (*Sn*) controlling for *M*, the unstandardised mean age at childbearing. These regression models are used for estimating adult mortality (Timæus 2013).

2. the relationship in single-year age groups between parental survival and summary indices of adult survivorship (*Lm*+*n*/*lm*) given the timing of childbearing (*M*). These models are useful when estimating fertility using the Own-Children Method as the estimate of orphans can be removed from the count of children of each age living apart from their mother, bringing the numerators and denominators of the fertility rates into correspondence.

The Devilish Stata do file applies both the conventional Own-Children Method and this enhanced version of the method to an IPUMS census file and compares the results (Timæus 2021).

By selecting the set of stable populations that Timæus (1992) adopted to produce coefficients for estimating adult mortality by the orphanhood method, the R script can be used to replicate the published estimation equations first calculated using a spreadsheet (*Supercalc*) in 1989.

The R script can also be used to calculate a set of multipliers for converting the proportions of mothers or fathers surviving to measures of life table survivorship or vice versa in a population with known characteristics (in particular, the age-pattern of fertility).

[![DOI](https://zenodo.org/badge/372783504.svg)](https://zenodo.org/badge/latestdoi/372783504)

#### References

Timæus, I. M. Estimation of adult mortality from paternal orphanhood: a reassessment and a new approach. *Population Bulletin of the United Nations*, 1992, **33**, 47–63. <https://blogs.lshtm.ac.uk/iantimaeus/files/2012/04/Timaeus-Pop-Bulletin-UN-Orphanhood-Paper.pdf>

Timæus, I. M. ‘Indirect estimation of adult mortality from orphanhood’, in *Tools for Demographic Estimation*, T. A. Moultrie, R. E. Dorrington, A. G. Hill, K. Hill, I. M. Timæus and B. Zaba (eds). International Union for the Scientific Study of Population, pp. 222–243, 2013. <https://demographicestimation.iussp.org/content/orphanhood>

Timæus, I. M. The Own-Children Method of fertility estimation - The devil is in the detail. *Demographic Research*, 2021, **45**(25), 825-840. DOI: [10.4054/DemRes.2021.45.25](<https://doi.org/10.4054/DemRes.2021.45.25>)
