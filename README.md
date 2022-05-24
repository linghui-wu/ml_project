# ECON31703 Topics in Econometrics Project

### Spatially Smooth Fixed Effects with Fused PPML

Reigner Kane, Linghui Wu

### Repository Structure

* `Code/` contains source code of the project

	* `*.R` constructs the adjacency matrics for residence and workplace census tracts, modifies the `pglasso` and `genlasso` packages, simulated the Monte Carlo DGP, estimates the fused lasso model and visualizes the estimates.
	
	* `*.jl` and `*.sbatch` estimates the fused ridge model.
	
	* `*.do` performs PPML gravity model estimation and draws maps of estimated fixed effects using the observed 2010 NYC commuting data.

* `Data/` 
	
	* `simulated_data` stores data from Monte Carlo simulations where fixed effects in the DGP are constant and independently drawn from a normal distribution.
	
	* `nyc2010_lodes_data` used for real world application exceeds the GitHub limitation for the file size and is available upon request.

* `Output/` contains point estimates from the fused PPML models with optimal penalty term.

* `Plots/` 
	
	* `RealData` includes NYC maps of the estimated origin and destination fixed effects from real data.
	
	* `SimulatedData` includes figures that visualize the origin and destination fixed effects from the generated data.

