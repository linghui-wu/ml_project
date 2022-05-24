clear all 

// Compute ppml FE estimates for full sample
use "../Data/nyc2010_lodes_wzero_wdelta.dta", clear
ppmlhdfe X_ij log_delta, absorb(fe_i_ppml=i fe_j_ppml=j) vce(cluster i j) d 
label var fe_i_ppml "Residence tract fixed effects"
label var fe_j_ppml "Workplace tract fixed effects"
tempfile df
save `df', replace
collapse (firstnm) fe_j_ppml, by(j)
export delimited "../Data/ppml_fe_j.csv", replace 
use i fe_i_ppml using `df', clear
collapse (firstnm) fe_i_ppml, by(i)
export delimited "../Data/ppml_fe_i.csv", replace

// Compute ppml FE estimates for sumsample (Bronx and Manhattan)
use `df', clear
gen county_i = substr(i, 1, 5)
gen county_j = substr(j, 1, 5)
keep if county_i == "36061" || county_i == "36085"
keep if county_j == "36061" || county_j == "36085"
drop fe_*
ppmlhdfe X_ij log_delta, absorb(fe_i_ppml=i fe_j_ppml=j) vce(cluster i j) d 
label var fe_i_ppml "Residence tract fixed effects"
label var fe_j_ppml "Workplace tract fixed effects"
tempfile df_subsample
save `df_subsample', replace 
