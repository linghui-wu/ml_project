clear all

graph set window fontface "Garamond"
graph set eps fontface "Times" 

cap program drop map_cutoff
cap program define map_cutoff
syntax, var(string) format(string) ///
    var_type(string) pct_list(numlist) output(string) fc(string)


// Calculate the percentiles for each observation
bysort `var': egen count_`var' = count(`var') 
drop if missing(count_`var')
gen cum_count_`var' = count_`var'[1]
replace cum_count_`var' = count_`var'[_n] + cum_count_`var'[_n-1] if _n > 1
gen pct_`var' = cum_count_`var' * 100 / cum_count_`var'[_N]

// Obtain sum stats for non-treated tracts
sum `var'
local min_`var' = r(min)
local max_`var' = r(max)

// Find the corresponding values for percentiles in `pct_list' and generate labels
local list_`var' = ""
local i 0 // Used to store the length of `pct_list'

foreach pct in `pct_list' {
    local i = `i' + 1

    local p`i' = string(`pct', "`format'")
    local pp`i' = "p`p`i''"

    qui count if pct_`var' < `pct'
    local val`i' = `var'[`r(N)']
    local list_`var' = "`list_`var''" + " " + "`val`i''"
}

local legend_list = ""
foreach ind of numlist 1/`i' {
    local ind_right = `ind' + 1

    if (`ind' == 1) {
        local min_`var' = string(`min_`var'', "`format'")
        local val`ind' = string(`val`ind'', "`format'")
        local temp_legend = "lab(`ind_right' min(`min_`var'') - " + ///
                            "`pp`ind''" + "(`val`ind''))"
    } 
    else {
        local ind_left = `ind' - 1
        local val`ind_left' = string(`val`ind_left'', "`format'")
        local val`ind' = string(`val`ind'', "`format'")
        local temp_legend = "lab(`ind_right' "+ "`pp`ind_left''" + ///
                            " - " +  "`pp`ind''" + "(`val`ind''))"
    }
    local legend_list = "`legend_list'" + " " + "`temp_legend'"

    if (`ind' == `i' ) {
        local ind_right = `ind' + 2
        local val`ind' = string(`val`ind'', "`format'")
        local max_`var' = string(`max_`var'', "`format'")
        local temp_legend_max = "lab(`ind_right' " +  "`pp`ind''" + ///
                            "(`val`ind'') - max" + "(`max_`var''))"
        local legend_list = "`legend_list'" + " " + "`temp_legend_max'"
    } 
}

maptile `var', geo(geoid11) geofolder("../Data/US_tract_2010_nyc") ///
    fc(`fc') res(0.3) ///
    cutv(`"`list_`var''"') ///
    twopt(legend(size(small) ///
        lab(1 "No `var_type' in 2010") ///
        `legend_list' )) ///
    savegraph(`output') replace 

end

// Plot maps for baseline estimates
import delimited "../Data/ppml_fe_i.csv", clear
rename i geoid11
local percentile_list 5 10 25 50 75 90 95
map_cutoff, var(fe_i_ppml) format("%10.1fc") ///
    var_type("estimates origin FEs") pct_list(`percentile_list') ///
    output("../Output/map_ppml_fe_i.png") fc(Purples)

import delimited "../Data/ppml_fe_j.csv", clear
rename j geoid11
local percentile_list 5 10 25 50 75 90 95
map_cutoff, var(fe_j_ppml) format("%10.1fc") ///
    var_type("estimates destination FEs") pct_list(`percentile_list') ///
    output("../Output/map_ppml_fe_j.png") fc(Purples)


// Plot maps for fused ridge estimates
import delimited "../Data/Estimated FEs/full_100_i.csv", clear
rename i geoid11
local percentile_list 5 10 25 50 75 90 95 
map_cutoff, var(fe) format("%10.1fc") ///
    var_type("estimated origin FEs") pct_list(`percentile_list') ///
    output("../Output/map_lasso_fe_i.png") fc(Purples)

import delimited "../Data/Estimated FEs/full_100_j.csv", clear
rename j geoid11
local percentile_list 10 25 50 75 90 95
map_cutoff, var(fe) format("%10.1fc") ///
    var_type("estimated destination FEs") pct_list(`percentile_list') ///
    output("../Output/map_lasso_fe_j.png") fc(Purples)


// Plot maps for fused ridge estimates
import delimited "../Output/fused_ridge_fe_i_optimal.csv", clear
rename i geoid11
local percentile_list 5 10 25 50 75 90 95
map_cutoff, var(fe_i_ridge_optimal) format("%10.1fc") ///
    var_type("estimated origin FEs") pct_list(`percentile_list') ///
    output("../Output/map_ridge_fe_i.png") fc(Purples)

import delimited "../Output/fused_ridge_fe_j_optimal.csv", clear
rename j geoid11
local percentile_list 5 10 25 50 75 90 95
map_cutoff, var(fe_j_ridge_optimal) format("%10.1fc") ///
    var_type("estimated destination FEs") pct_list(`percentile_list') ///
    output("../Output/map_ridge_fe_j.png") fc(Purples)
