/*
** Estimation of age-specific fertility by the own-children method **

This script applies the own-children method of estimating fertility to a pair
of successive censuses conducted in the same country and compares the results.
It computes cohort-period fertility rates for single years of age and 
five-year periods of time. Data on children aged 0 and 1 are discarded.

Two estimators are used - the conventional one proposed by Cho et al. (1986)
and a new estimator that brings the numerators and denominators of the rates
into correspondence by removing an estimate of orphaned children from the 
absentee (unmatched) children in the numerators rather than by adding an 
estimate of the women who have died since giving birth to the denominators.

The script is designed for use with IPUMS microdata files in which children 
have already been linked to their coresident mothers.

The final plot is used as Figure 2 of "The Own-Children Method of fertility 
estimation – the devil is in the detail".

Ian Timaeus
 
*! version 1.0 07-Aug-2021
***************************************************************************** */
version 13.0
set more off
clear all
local min_age = 15
local max_age = 79
cd "C:/Users/ecpsitim/temp/ipums"
/* Parameters for the USA 1900 and 1910
local dpath "C:/Users/ecpsitim/Documents/Ian - static/Old DataBases/Stata/"
local originalfile "usa_00002.dta"
global cc "US"
local censuses "1900 1910"
* Intercept and slope of alpha defining relational logit model life tables for
* intercensal period & period before 1st census - children then adults
matrix cm_coeff = (-0.065, 0.0108, -0.186, 0.0121)
matrix am_coeff = (-0.065, 0.0108, -0.186, 0.0121)
*/
* Parameters for Zambia 2000 and 2010
local dpath "C:/Users/ecpsitim/temp/ipums/"
local originalfile "ipumsi_00001"
local country = 894
global cc "ZM"
local censuses "2000 2010"
* Intercept and slope of alpha for intercensal period & before 1st census
matrix cm_coeff = (-0.15, 0.055, 0.4, 0.02)
matrix am_coeff = (0.6, 0.1, 1.6, -0.07)
* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
*
* MATA FUNCTION
mata:
mata clear
// =============================================================================
void initialise_mortality(real scalar cen, real scalar cen2, real scalar mbar)
//	Set up matrices of estimates of child & adult mortality & maternal survival
{
	real matrix  lx, Lx, Px, ms_coeff
	real colvector  Ysx, dead, m_alive, gone
	real rowvector  cm_coeff, am_coeff, l25
	real scalar  topr, i, len

	cm_coeff = st_matrix("cm_coeff")
	am_coeff = st_matrix("am_coeff")
	topr = strtoreal(st_local("max_age")) - 11
//	Initialise a 1-parameter Princeton West logit relational model (e0=50, both 
//	sexes, SRAB 1.035) of survivorship by age (0->100)
	Ysx =   (-1.9862\ -1.7495\ -1.6396\ -1.5736\ -1.5288\ -1.4981\ -1.4748\
	-1.4558\ -1.4390\ -1.4232\ -1.4081\ -1.3935\ -1.3791\ -1.3643\ -1.3486\
	-1.3314\ -1.3128\ -1.2929\ -1.2718\ -1.2495\ -1.2261\ -1.2017\ -1.1769\
	-1.1518\ -1.1269\ -1.1022\ -1.0776\ -1.0529\ -1.0282\ -1.0033\ -0.9781\
	-0.9526\ -0.9269\ -0.9010\ -0.8748\ -0.8483\ -0.8214\ -0.7942\ -0.7666\
	-0.7385\ -0.7100\ -0.6809\ -0.6513\ -0.6213\ -0.5907\ -0.5598\ -0.5282\
	-0.4959\ -0.4624\ -0.4274\ -0.3909\ -0.3527\ -0.3130\ -0.2719\ -0.2294\
	-0.1857\ -0.1404\ -0.0932\ -0.0436\ 0.0088\ 0.0643\ 0.1228\ 0.1844\ 0.2488\
	0.3160\ 0.3860\ 0.4592\ 0.5360\ 0.6169\ 0.7025\ 0.7931\ 0.8890\ 0.9906\
	1.0984\ 1.2128\ 1.3342\ 1.4633\ 1.6005\ 1.7467\ 1.9025\ 2.0687\ 2.2462\
	2.4358\ 2.6384\ 2.8550\ 3.0864\ 3.3338\ 3.5980\ 3.8802\ 4.1813\ 4.5023\
	4.8441\ 5.2075\ 5.5933\ 6.0020\ 6.4342\ 6.8902\ 7.3702\ 7.8744\ 8.4028)
	len = rows(Ysx)
//	Invert life table - place the most recent births in the bottom row	
	Ysx = Ysx[len..1]

//	DEATHS OF CHILDREN: Produce cohort estimates of Q(x)

	lx = J(len, len, .)
	Lx = J(len, len, .)
	Px = J(len, len-1, .)
	for (i=1; i<=len; i++) lx[, i] = 1 :/ (J(len, 1, 1) + 
		exp(cm_coeff[1] :+ cm_coeff[2] :* min((20, i-1)) :+ Ysx[1..len]))
	if (cen==cen2) {
	    lx[, 11..len] = lx[, 1..len-10]
		for (i=1; i<=10; i++) lx[, i] = 1 :/ (J(len, 1, 1) + 
			exp(cm_coeff[3] :+ cm_coeff[4] :* (i-1) :+ Ysx[1..len]))
	}
//	Convert lx to Lx within each period life table for an exact date
	Lx[1..len-1, ] = (lx[1..len-1, ] + lx[2..len, ]) :/2
	Lx[len, ] = 0.3 :+ 0.7 :* lx[len, ]
//	Estimate single-year cohort survivorship over time
	Px[1..len-1, ] = ((Lx[1..len-1, 1..len-1] :/ Lx[2..len, 1..len-1]) + 
		(Lx[1..len-1, 2..len] :/ Lx[2..len, 2..len])) :/ 2
	Px[len, 1..len-1] = (Lx[len, 1..len-1] + Lx[len, 2..len]) / 2		
//	Work leftward calculating cohort Lx from the Px along the diagonals
	for (i=len-1; i>=1; i--) Lx[1..len-1, i] = Px[1..len-1, i] :* 
		Lx[2..len, i+1]
	Lx[len, 1..len-1] = Px[len, ]
//	Col. 1 now holds survivorship (Lx) for cohorts of children by age at census 
	dead = 1 :- Lx[len+1-topr..len, 1]
	st_matrix("dead", dead)
	
//	OLD METHOD - women's survivorship (same models as for children in the USA)

	for (i=1; i<=len; i++) lx[, i] = 1 :/ (J(len, 1, 1) + 
		exp(am_coeff[1] :+ am_coeff[2] :* min((20, i-1)) :+ Ysx[1..len]))
	if (cen==cen2) {
	    lx[, 11..len] = lx[, 1..len-10]
		for (i=1; i<=10; i++) lx[, i] = 1 :/ (J(len, 1, 1) + 
			exp(am_coeff[3] :+ am_coeff[4] :* (i-1) :+ Ysx[1..len]))
	}
//	Convert lx to Lx within each period life table for an exact date
	Lx[1..len-1, ] = (lx[1..len-1, ] + lx[2..len, ]) :/2
	Lx[len, ] = 0.3 :+ 0.7 :* lx[len, ]
//	Estimate single-year cohort survivorship over time
	Px[1..len-1, ] = ((Lx[1..len-1, 1..len-1] :/ Lx[2..len, 1..len-1]) + 
		(Lx[1..len-1, 2..len] :/ Lx[2..len, 2..len])) :/ 2
	Px[len, 1..len-1] = (Lx[len, 1..len-1] + Lx[len, 2..len]) / 2		
//	Work leftward calculating cohort Lx from the Px along the diagonals
	for (i=len-1; i>=1; i--) Lx[1..len-1, i] = Px[1..len-1, i] :* 
		Lx[2..len, i+1]
	Lx[len, 1..len-1] = Px[len, ]
//	Transfer to Stata for calculating mothers' survivorship in the "old" OCM
	st_matrix("Lx", Lx)

//	NEW METHOD - estimate the prevalence of orphanhood by age of child
	
	l25 = (Lx[75, ] + Lx[76, ]) :/ 2
//	Coefficients to predict proportion of children with living mothers
	ms_coeff =	(0.1242, -0.000059, 0.8773\ 0.0097, -0.000184, 0.9953\
				-0.0130, -0.000301, 1.0212\ -0.0253, -0.000416, 1.0365\
				-0.0366, -0.000532, 1.0509\ -0.0455, -0.000655, 1.0631\
				-0.0505, -0.000787, 1.0716\ -0.0538, -0.000931, 1.0787\
				-0.0555, -0.001089, 1.0846\ -0.0562, -0.001262, 1.0899\
				-0.0567, -0.001452, 1.0954\ -0.0559, -0.001659, 1.1000\
				-0.0538, -0.001885, 1.1039\ -0.0507, -0.002130, 1.1073\
				-0.0470, -0.002396, 1.1105\ -0.0426, -0.002684, 1.1136\
				-0.0377, -0.002994, 1.1167\ -0.0316, -0.003327, 1.1193\
				-0.0243, -0.003685, 1.1211\ -0.0161, -0.004069, 1.1228\
				-0.0064, -0.004481, 1.1236\ 0.0045, -0.004921, 1.1239\
				0.0163, -0.005391, 1.1238\ 0.0297, -0.005892, 1.1230\
				0.0445, -0.006425, 1.1213\ 0.0608, -0.006992, 1.1190\
				0.0785, -0.007593, 1.1159\ 0.0977, -0.008228, 1.1121\
				0.1185, -0.008900, 1.1072\ 0.1408, -0.009607, 1.1015\
				0.1647, -0.010349, 1.0948\ 0.1903, -0.011127, 1.0870\
				0.2174, -0.011939, 1.0782\ 0.2460, -0.012784, 1.0683\
				0.2759, -0.013658, 1.0573\ 0.3071, -0.014558, 1.0452\
				0.3393, -0.015478, 1.0321\ 0.3722, -0.016414, 1.0179\
				0.4058, -0.017357, 1.0027\ 0.4395, -0.018298, 0.9867\
				0.4731, -0.019227, 0.9698\ 0.5060, -0.020132, 0.9523\
				0.5379, -0.021000, 0.9343\ 0.5683, -0.021818, 0.9159\
				0.5966, -0.022571, 0.8974\ 0.6223, -0.023244, 0.8789\
				0.6451, -0.023823, 0.8607\ 0.6644, -0.024292, 0.8429\
				0.6798, -0.024638, 0.8257\ 0.6909, -0.024849, 0.8095\
				0.6975, -0.024915, 0.7943\ 0.6993, -0.024826, 0.7805\
				0.6962, -0.024576, 0.7682\ 0.6880, -0.024162, 0.7574\
				0.6748, -0.023583, 0.7485\ 0.6569, -0.022843, 0.7415\
				0.6342, -0.021948, 0.7364\ 0.6072, -0.020912, 0.7334\
				0.5763, -0.019748, 0.7327\ 0.5419, -0.018478, 0.7346\
				0.5046, -0.017126, 0.7396\ 0.4651, -0.015714, 0.7499\
				0.4241, -0.014272, 0.7667\ 0.3822, -0.012828, 0.7917\
				0.3404, -0.011405, 0.8288\ 0.2993, -0.010022, 0.8852\
				0.2597, -0.008702, 0.9668\ 0.2222, -0.007458, 1.0839)
// 	Estimate the proportions by age with living mothers
	m_alive = ms_coeff[topr..1, 1] + (ms_coeff[topr..1, 2] :* 
		mbar) + (ms_coeff[topr..1, 3] :* Lx[76-topr..75, 1] :/ l25[topr..1]')
//	Probability that an unorphaned child is not with their mother
	gone = st_matrix("prop_gone")
	gone = gone[, 3]
	gone = 1 :- (1 :- gone) :/ m_alive[1..topr]
	gone = (gone :> J(topr,1, 0)) :* gone
	gone = (gone :<= J(topr,1, 1)) :* gone :+ (gone :> J(topr,1, 1))
//	prop_gone is now estimated absentee children with living mothers 
	st_matrix("adj_prop_gone", gone)
}
end
* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
*
* Set up file of women aged 15-79
*
use "`dpath'`originalfile'", clear
capture drop if country~=`country'
keep year serial pernum perwt sex age chborn chsurv
keep if sex==2 & inrange(age, `min_age',`max_age')
compress
gen double mth_id = serial * 100 + pernum
* Normalise weights on women
tempvar meanwt
bys year: egen `meanwt' = mean(perwt)
gen wt = perwt/`meanwt'
replace wt = 1 if wt==.
* The data are classified by cohort of women defined by age at census
rename age coh
* Rename summary birth history variables
gen byte tceb = chborn - ("$cc"=="US") if chborn>0
gen byte tcs = chsurv - ("$cc"=="US") if chsurv>0
gen byte tcd = tceb - tcs
sort year mth_id
save ${cc}_womfile, replace
*
* Set up file of own children
*
use "`dpath'`originalfile'", clear
capture drop if country~=`country'
keep year serial age momloc stepmom
keep if momloc>0 & stepmom==0
compress
rename age intage
* Add variables from mother's record to their children's records
gen double mth_id = serial * 100 + momloc
drop serial momloc stepmom
sort year mth_id intage
merge m:1 year mth_id using ${cc}_womfile, nolabel
keep if _merge==3
sort year mth_id intage
save ${cc}_birthfile, replace
*
* Calculate smoothed proportions elsewhere
*
use "`dpath'`originalfile'", clear
capture drop if country~=`country'
keep if age<(`max_age' - 11)
gen byte non_own = momloc==0 | inrange(stepmom, 1,7)
gen byte pop = 1
gcollapse (sum) non_own pop [pw=perwt], by(year age)
gen paway = .
* Smooth curve by means of logistic regression with cubic splines
foreach year of local censuses {
	mkspline spa = age if year==`year', cubic nknots(5)
	glm non_own spa* if year==`year', family(binomial pop)
	predict p_away if year==`year'
	replace paway = p_away/pop if year==`year'
	drop spa* p_away
}
gsort year -age
keep year age paway
save _${cc}_prop_kids_away, replace
*
* Estimate child mortality (in this version, these estimates guide the selection 
* of model life tables rather than get directly incorporated in the analysis)
*
foreach year of local censuses {
	use ${cc}_womfile.dta, clear
	keep if year==`year'
	gen byte idx = floor(coh/5)-2
	* Drop women aged <15 or 70+
	drop if ~inrange(idx, 1,11)
	gen byte nwom = 1
	* Collapse to file on 5-year age groups by region & level of schooling
	drop __00*
	gcollapse (sum) tceb tcd nwom [pw=wt], by(idx)
	local p1p2 = (tceb[_n]/nwom[_n])/(tceb[_n+1]/nwom[_n+1])
	local p2p3 = (tceb[_n+1]/nwom[_n+1])/(tceb[_n+2]/nwom[_n+2])
	* Coefficients for multipliers to convert proportions dead to q(x)
	* (West coefficients, data tabulated by mother's age group)
	matrix coeffs = (1.1415, -2.7070,  0.7663 \ ///
					 1.2563, -0.5381, -0.2637 \ ///
					 1.1851,  0.0633, -0.4177 \ ///
					 1.1720,  0.2341, -0.4272 \ ///
					 1.1865,  0.3080, -0.4452 \ ///
					 1.1746,  0.3314, -0.4537 \ ///
					 1.1639,  0.3190, -0.4435)
	* Coefficients for estimating the time location of mortality measures
	matrix tcoeffs =(1.097,   5.5628, -1.9956 \ ///
					1.3062,   5.5677,  0.2962 \ ///
					1.5305,   2.5528,  4.8962 \ ///
					1.9991,  -2.4261, 10.4282 \ ///
					2.7632,  -8.4065, 16.1787 \ ///
					4.3468, -13.2436, 20.199  \ ///
					7.5242, -14.2013, 20.0162)
	* Standard life table for children (West, e0=50, sexes combined, SRAB 1.035)
	* Y(2), Y(3), Y(5), Y(10) ... Y(55)
	matrix Ysx = (-0.87474\ -0.81982\ -0.76440\ -0.71160\ -0.67428\ ///
		 -0.62476\ -0.56344\ -0.50163\ -0.43739\ -0.36927\ -0.29537\ ///
		 -0.21372\ -0.11472)
	* Attach appropriate logit from standard to records for cohorts of mothers
	gen Ys = 2 * Ysx[idx,1]
	* Index to remove children of teenagers from computation of mortality trend
	gen byte teen = idx==1
	*
	* Overall estimates of mortality relative to the standard (2*alpha)
	*
	* Divide births by Brass multipliers (equivalent to adjusting proportions  
	* dead to q(x) but leaves counts of events unchanged for se's)
	gen adj_tceb = tceb / (coeffs[idx, 1] + coeffs[idx, 2]*`p1p2'  ///
		+ coeffs[idx, 3]*`p2p3') if idx<=7
	* Approximate multipliers for woman>50 based on original Brass multipliers
	qui replace adj_tceb = tceb / (0.15 + coeffs[7,1] + coeffs[7,2]*`p1p2' ///
		+ coeffs[7,3]*`p2p3') if idx>7
	replace adj_tceb = tcd + 0.0001 if tcd>=adj_tceb
	* Estimate time locations of the mortality estimates
	gen tloc = tcoeffs[idx,1] + tcoeffs[idx,2]*`p1p2'  ///
		+ tcoeffs[idx,3]*`p2p3' if idx<=7
	* Once childbearing is over, children's exposure rises with mother's age
	qui replace tloc = tcoeffs[7,1] + tcoeffs[7,2]*`p1p2'  ///
		+ tcoeffs[7,3]*`p2p3' + (idx-7)*5 if idx>7
	* Time-located Brass estimates of 2*alpha relative to the standard
	glm tcd c.tloc i.teen if adj_tceb>0 & inrange(idx, 1,11) [aw=nwom],  ///
		offset(Ys) fam(binomial adj_tceb) link(logit)
	matrix tmort_ests = get(_b)
	matrix mt_coeff = (tmort_ests[1,4], tmort_ests[1,1])
	clear
	svmat mt_coeff
	save _${cc}_mt_coeff`year', replace
}
*
* Calculate matrices of cohort-period fertility rates
*
tokenize `censuses'
capture erase "${cc}_OCM.dta"
foreach year of local censuses {
	* Calculate mean age at childbearing from own-children aged 0
	use "${cc}_birthfile" if intage==0 & year==`year', clear
	su coh, meanonly
	local mbar = r(mean) - 0.5
	* Transfer proportions of children elsewhere & who have died to matrices
	use "_${cc}_prop_kids_away" if year==`year', clear
	mkmat * , matrix(prop_gone)
	mata:initialise_mortality(`year', `2', `mbar')
	use "${cc}_birthfile" if year==`year', clear
	keep mth_id coh intage wt
	gsort mth_id -intage
	* Assign children to their mean age, spacing multiple births by a day 
	gen age = intage + 0.5 - ((intage==intage[_n-1] & mth_id==mth_id[_n-1]) ///
		+ (intage==intage[_n-2] & mth_id==mth_id[_n-2]))/365	
	gen byte closed = 1
	* Create records for each parous woman's open interval from the woman's file
	tempfile temp_ocm
	save `temp_ocm', replace
	use mth_id tceb coh wt year using "${cc}_womfile" if year==`year' , clear
	drop if tceb==0
	drop tceb
	gen age = 0
	gen byte closed = 0
	append using `temp_ocm'
	gsort mth_id -closed -age
	* Calculate mother's age at the birth of each of her children
	gen m_ageb = coh + 0.5 - age
	* stset as a repeatable event indexed by age at birth (for ASFRs)
	stset m_ageb, failure(closed) id(mth_id) exit(time .)
	* Drop childhood & post-menopausal exposure to reduce file size
	qui stsplit edges, at(11.5 49.5)
	replace closed = 0 if closed==. & edges[_n+1]==49.5 & mth_id==mth_id[_n+1]
	replace age = . if ~closed
	qui drop if edges ~= 11.5
	qui drop edges
	* Split exposure into yearly *periods* before census with dummies
	stsplit i_per, every(1) after(time = 11.5)
	replace i_per = coh - 12 - i_per
	gen exposure = _t - _t0
	* Collapse data into APC table
	capture drop __0000*
	gcollapse (sum) b=_d exposure [pw=wt], by(i_per coh)
	gen adj_b = b / (1 - dead[`max_age' - 11 - i_per, 1])
	* Old OCM estimator
	gen adj_b_old = adj_b / (1 - prop_gone[`max_age' - 11 - i_per, 3]) * ///
		2 * Lx[100-coh, 1] / (Lx[101-coh+i_per, i_per+1] ///
		+ Lx[100-coh+i_per, i_per+1])
	* The devil is in the detail	
	gen adj_b_new = adj_b / (1 - adj_prop_gone[`max_age' - 11 - i_per, 1])
	gen int census_year = `year'
	* Recreate AGE in collapse COHORT-PERIOD file
	gen byte m_ageb = coh - i_per
	* Indicators for 5-year periods discarding births in the last 2 years
	gen int per5 = census_year - (7+5*floor((i_per-2)/5))
	foreach t of numlist 62 (-5) 7 {
		local startyr = census_year - `t'
		local endyr = census_year - `t' + 5
		label def fiveyr `startyr' "`startyr'–`endyr'", modify
	}
	label val per5 fiveyr
	capture append using "${cc}_OCM.dta"
	save "${cc}_OCM.dta", replace
}
*
* Plot TFRS for each method and census
*
* Collapse into rates for single-year periods
collapse (sum) b adj_b adj_b_new adj_b_old exposure, ///
	by(census_year i_per m_ageb)
gen rate_new = adj_b_new/exposure
gen rate_old = adj_b_old/exposure
replace i_per = census_year-i_per
* Calculate single-year period TFRs from the ASFR
bys census_year i_per: egen tfr_n = total(rate_new)
bys census_year i_per: egen tfr_o = total(rate_old)
tokenize `censuses'
tw (line tfr_n i_per if census_year==`1', ///
		lcolor(erose) lpattern(shortdash)) ///
	(line tfr_n i_per if census_year==`2', ///
		lcolor(eltgreen) lpattern(shortdash)) ///
	(line tfr_o i_per if census_year==`1', color(erose)) ///
	(line tfr_o i_per if census_year==`2', color(eltgreen)) ///
	if inrange(census_year-i_per, 7,22), plotregion(lstyle(none)) ///
	ytitle(Total fertility, xoffset(-2)) legend(rows(2) nobox  ///
	region(lstyle(none)) label(1 "`1' New") label(2 "`2' New") ///
	label(3 "`1' Old") label(4 "`2' Old") label(5 "`1' Mor") ///
	label(6 "`2' Mor") label(7 "`1' Own") label(8 "`2' Own"))
*
* Plot ASFRs for each method and census (Figure 2)
*
use "${cc}_OCM.dta", clear
* Collapse into rates for 5 year periods
collapse (sum) b adj_b adj_b_new adj_b_old exposure, by(census_year per5 m_ageb)
gen rate_new = adj_b_new/exposure
gen rate_old = adj_b_old/exposure
bys census_year per5: egen tfr_n = total(rate_new)
bys census_year per5: egen tfr_o = total(rate_old)
* Prorate the overlapping rates from the 2nd census to the TFR in the 1st census
preserve
	keep if (census_year~=census_year[_n-1] | per5~=per5[_n-1]) & ///
		inrange(census_year - per5, 7, 22)	
	gsort -per5 census_year
	mkmat per5 census_year tfr_n tfr_o, matrix(tfrs)
restore
egen int date1 = min(census_year)
egen int date2 = max(census_year)
* Only adjust the rates estimated for 12-22 years before the later census
replace rate_new = rate_new * tfrs[3,3]/tfrs[4,3] if per5==(date1 - 7) ///
	& census_year~=date1
replace rate_new = rate_new * tfrs[5,3]/tfrs[6,3] if per5==(date1 - 12) ///
	& census_year~=date1
replace rate_old = rate_old * tfrs[3,4]/tfrs[4,4] if per5==(date1 - 7) ///
	& census_year~=date1
replace rate_old = rate_old * tfrs[5,4]/tfrs[6,4] if per5==(date1 - 12) ///
	& census_year~=date1
local d1 = date1[1]
local d2 = date2[1]
tw (line rate_new m_ageb if census_year==date1, lcolor(cranberry) ///
		lpatt(shortdash)) ///
	(line rate_old m_ageb if census_year==date1, lcolor(erose)) ///
	(line rate_new m_ageb if census_y~=date1, lcolor(midblue) ///
		lpatt(shortdash)) ///
	(line rate_old m_ageb if census_year~=date1, lcolor(eltgreen)) ///
	if inrange(census_year-per5, 7,22), plotregion(lstyle(none)) ///
 	ytitle(Fertility rate)  yscale(range(0 0.3)) ylabel(, labsize(medium)) ///
	by(per5, imargin(small) note("") cols(3) b1title(Age, size(medsmall))) ///
	legend(rows(1) colfirst nobox region(lstyle(none)) label(1 "`d1' – New") ///
	label(2 "`d1' – Old") label(3 "`d2' – New") label(4 "`d2' – Old")) ///
	xlabel(, labsize(medium))


