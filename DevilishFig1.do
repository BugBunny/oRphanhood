/*
** Figure 1 of "The Own-Children Method of fertility estimation – the devil
** is in the detail"

Plots the proportion of the children reported in the birth histories collected 
in 9 DHS surveys that are living apart from their mother over her age group by
the children's age group.

Ian Timaeus
*! version 1.0 07-Aug-2021
*******************************************************************************/
version 13.0
global mydir "C:\Users\ecpsitim\Documents\Ian - static\000DHS FL files\"
cd "$mydir"
tempfile temp1
local cfiles "AZIR52FL EGIR51FL CMIR44FL HTIR52FL IAIR52FL"
local cfiles "`cfiles' PEIR51FL PHIR52FL SNIR4HFL TZIR4QFL"

foreach cntry of local cfiles {
	use caseid v005 v013 bidx* b5* b8* b16* using `cntry', clear
	* v005 - Sample weight
	* v013 - Five-year age group
	* b5 - Child is alive?
	* b8 - Child's age
	* b16 - Child's line number in household member file
	rename b*_0* b*_*
	reshape long bidx_ b5_ b8_ b16_, i(caseid)
	drop if bidx_==. | b8>=20
	rename *_ *
	gen byte non_own = b16==0 if b5==1
	gen wt = v005/1000000
	replace b8 = 5*int(b8/5)
	replace v013 = v013*5+12.5
	gen nwom = 1/wt
	collapse (mean) non_own (sum) nwom [pw=wt], by(v013 b8)
	drop if nwom<25
	gen country = substr("`cntry'",1,2)	
	capture append using `temp1'
	save `temp1', replace
}
replace country = "Azerbaijan, 2006" if country=="AZ"
replace country = "Cameroon, 2004" if country=="CM"
replace country = "Egypt, 2005" if country=="EG"
replace country = "Haiti, 2005" if country=="HT"
replace country = "India, 2005" if country=="IA"
replace country = "Peru, 2004" if country=="PE"
replace country = "Philippines, 2008" if country=="PH"
replace country = "Senegal, 2005" if country=="SN"
replace country = "Tanzania, 2004" if country=="TZ"
twoway (line non_own v013 if inrange(b8,0,4), sort) (line non_own v013 ///
	if inrange(b8,5,9), sort) (line non_own v013 if inrange(b8,10,14), sort) ///
	(line non_own v013 if inrange(b8,15,19), sort), legend(rows(1) ///
	title("Age of child", size(small)) size(small) region(lstyle("none"))	///
	label(1 "0–4") label(2 "5–9") label(3 "10–14") label(4 "15–19")) ///
	ytitle("Proportion of children not with mother", size(small)) ///
	xtitle("Age of mother", size(small)) plotregion(lstyle("none")) ///
	by(country, note("")) xlabel(15(10)45, labsize(3)) ///
	ylabel(0(.1).5, labsize(3)) subtitle(,size(small)) 
