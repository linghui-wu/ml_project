using Base.Threads, CSV, DataFrames, JLD, Optim

gamma = parse(Float64,ARGS[1])

@time data = sort(DataFrame(load("../Data/nyc2010_lodes_wzero_wdelta.dta")), [:j, :i])
num_i = length(unique(data.i))
num_j = length(unique(data.j))
i_df = DataFrame(id_i=collect(1:num_i), i=sort(unique(data.i)))
j_df = DataFrame(id_j=collect(1:num_j), j=sort(unique(data.j)))

fe_i_ppml = CSV.read("../Data/ppml_fe_i.csv", DataFrame).fe_i_ppml
fe_j_ppml = CSV.read("../Data/ppml_fe_j.csv", DataFrame).fe_j_ppml

neighboring_tracts = CSV.read("../Data/nlist_2010.csv", DataFrame)
rename!(neighboring_tracts, :SOURCE_TRACTID=>:i)
neighboring_tracts.i = string.(neighboring_tracts.i)
neighboring_tracts.NEIGHBOR_TRACTID = string.(neighboring_tracts.NEIGHBOR_TRACTID)
i_neighboring = innerjoin(i_df, neighboring_tracts, on=:i)
rename!(i_df, :i=>:NEIGHBOR_TRACTID)
i_neighboring = innerjoin(i_neighboring, i_df, on=:NEIGHBOR_TRACTID, makeunique=true)

neighboring_tracts = CSV.read("../Data/nlist_2010.csv", DataFrame)
rename!(neighboring_tracts, :SOURCE_TRACTID=>:j)
neighboring_tracts.j = string.(neighboring_tracts.j)
neighboring_tracts.NEIGHBOR_TRACTID = string.(neighboring_tracts.NEIGHBOR_TRACTID)
j_neighboring = innerjoin(j_df, neighboring_tracts, on=:j)
rename!(j_df, :j=>:NEIGHBOR_TRACTID)
j_neighboring = innerjoin(j_neighboring, j_df, on=:NEIGHBOR_TRACTID, makeunique=true)

ϵ = 7.98
params_ppml = vcat(ϵ, fe_i_ppml, fe_j_ppml)

function ll(params; γ=gamma, data=data, num_i=num_i, num_j=num_j, D_i=i_neighboring, D_j=j_neighboring)
   ϵ = params[1]
   α_i_vec = params[2:num_i + 1]
   α_j_vec = params[num_i+2:end]
   i_df = DataFrame(id_i=collect(1:num_i), i=unique(data.i));
   i_df.fe_i = α_i_vec; 
   j_df = DataFrame(id_j=collect(1:num_j), j=unique(data.j));
   j_df.fe_j = α_j_vec; 
   data = outerjoin(data, i_df, on=:i);
   data = outerjoin(data, j_df, on=:j);
   term1 = -data.X_ij .* (-ϵ .* data.log_delta .+ data.fe_i .+ data.fe_j);
   term2 = data.delta .^ -ϵ .* exp.(data.fe_i .+ data.fe_j);
   i_nbhd = outerjoin(i_df, D_i, on=:i, makeunique=true);
   rename!(i_nbhd, :i=>:i_origin, :fe_i=>:fe_i_origin, :NEIGHBOR_TRACTID=>:i)
   i_nbhd = outerjoin(i_df, i_nbhd, on=:i, makeunique=true, matchmissing=:equal)
   i_nbhd = i_nbhd[completecases(i_nbhd), :]  # Drop row with missing values in columns
   j_nbhd = outerjoin(j_df, D_j, on=:j, makeunique=true)
   rename!(j_nbhd, :j=>:j_origin, :fe_j=>:fe_j_origin, :NEIGHBOR_TRACTID=>:j)
   j_nbhd = outerjoin(j_df, j_nbhd, on=:j, makeunique=true, matchmissing=:equal)
   j_nbhd = j_nbhd[completecases(j_nbhd), :]  # Drop row with missing values in columns
   penalty_i = sum((i_nbhd.fe_i .- i_nbhd.fe_i_origin) .^ 2)
   penalty_j = sum((j_nbhd.fe_j .- j_nbhd.fe_j_origin) .^ 2)
   log_likelihood_val = sum(skipmissing(term1) .+ skipmissing(term2)) + γ * (penalty_i + penalty_j)
   println(ϵ, " ", log_likelihood_val)
    return log_likelihood_val
end

function g!(G, params; data=data, γ=gamma, num_i=num_i, num_j=num_j, D_i=i_neighboring, D_j=j_neighboring)
    ϵ = params[1]
    α_i_vec = params[2:num_i + 1]
    α_j_vec = params[num_i+2:end]
    i_df = DataFrame(id_i=collect(1:num_i), i=unique(data.i));
    i_df.fe_i = α_i_vec;
    j_df = DataFrame(id_j=collect(1:num_j), j=unique(data.j));
    j_df.fe_j = α_j_vec;
    data = outerjoin(data, i_df, on=:i);
    data = outerjoin(data, j_df, on=:j);
    i_nbhd = outerjoin(i_df, D_i, on=:i, makeunique=true)
    rename!(i_nbhd, :i=>:i_origin, :fe_i=>:fe_i_origin, :NEIGHBOR_TRACTID=>:i)
    select!(i_nbhd, Not(:id_i_2))
    i_nbhd = outerjoin(i_df, i_nbhd, on=:i, makeunique=true, matchmissing=:equal)
    # i_nbhd.diff_term = 2 * γ * (i_nbhd.fe_i_origin .- i_nbhd.fe_i)
    # i_nbhd = i_nbhd[completecases(i_nbhd), :]  # Drop row with missing values in columns
    j_nbhd = outerjoin(j_df, D_j, on=:j, makeunique=true)
    rename!(j_nbhd, :j=>:j_origin, :fe_j=>:fe_j_origin, :NEIGHBOR_TRACTID=>:j)
    select!(j_nbhd, Not(:id_j_2))
    j_nbhd = outerjoin(j_df, j_nbhd, on=:j, makeunique=true, matchmissing=:equal)
    # j_nbhd.diff_term = 2 * γ * (j_nbhd.fe_j_origin .- j_nbhd.fe_j)
    j_nbhd = j_nbhd[completecases(j_nbhd), :]  # Drop row with missing values in columns
    G[1] = sum(skipmissing(data.log_delta .* (data.X_ij .- exp.(data.fe_i .+ data.fe_j) .* data.delta .^ -ϵ)))
    Threads.@threads for idx in 1:num_i
        data_i = data[data.id_i .== idx, :];
        term1_i = sum(skipmissing(-data_i.X_ij + exp.(data_i.fe_i .+ data_i.fe_j) .* data_i.delta .^ -ϵ))
        i_nbhd_i = i_nbhd[.!ismissing.(i_nbhd.i_origin), :]
        i_nbhd_i = i_nbhd_i[i_nbhd_i.id_i_2 .== idx, :]
        term2_i = sum(2 * 2 * γ * (i_nbhd_i.fe_i_origin .- i_nbhd_i.fe_i))
        G[1 + idx] = sum(skipmissing(vcat(term1_i, term2_i)))
    end
    Threads.@threads for idx in 1:num_j
        data_j = data[data.id_j .== idx, :];
        term1_j = sum(skipmissing(-data_j.X_ij + exp.(data_j.fe_i .+ data_j.fe_j) .* data_j.delta .^ -ϵ))
        j_nbhd_j = j_nbhd[.!ismissing.(j_nbhd.j_origin), :]
        j_nbhd_j = j_nbhd_j[j_nbhd_j.id_j_2 .== idx, :]
        term2_j = sum(2 * 2 * γ * (j_nbhd_j.fe_j_origin .- j_nbhd_j.fe_j))
        G[1 + num_i + idx] = sum(skipmissing(vcat(term1_j, term2_j)))
    end
    return G
end


# Manual gradient check
# params1 = vcat(ϵ - 1e-3, fe_i_ppml, fe_j_ppml);
# params2 = vcat(ϵ + 1e-3, fe_i_ppml, fe_j_ppml);
# G1 = (ll(params2) - ll(params1)) / 2e-3

# params3 = vcat(ϵ, vcat(fe_i_ppml[1] - 1e-3, fe_i_ppml[2:end]), fe_j_ppml);
# params4 = vcat(ϵ, vcat(fe_i_ppml[1] + 1e-3, fe_i_ppml[2:end]), fe_j_ppml);
# G2 = (ll(params4) - ll(params3)) / 2e-3

# params5 = vcat(ϵ, fe_i_ppml, vcat(fe_j_ppml[1] - 1e-3, fe_j_ppml[2:end]));
# params6 = vcat(ϵ, fe_i_ppml, vcat(fe_j_ppml[1] + 1e-3, fe_j_ppml[2:end]));
# G2162 = (ll(params6) - ll(params5)) / 2e-3

# G = zeros(num_i + num_j + 1);
# @time g!(G, params_ppml);

# (G1, G[1])
# (G2, G[2])
# (G2162, G[2162])


# Optimization
lower = -10 * ones(num_i + num_j + 1); lower[1] = 6.0;
upper = 10 * ones(num_i + num_j + 1); upper[1] = 9.0;
@time results = optimize(ll, g!, lower, upper, params_ppml, Fminbox(GradientDescent()), 
        Optim.Options(show_trace=true, 
            x_abstol=1e-2, x_reltol=1e-2, f_tol=1e-2, g_tol=1e-2, iterations=1, 
            outer_x_abstol=1e-2, outer_f_abstol=1e-2, outer_g_abstol=1e-3, outer_iterations=1))
estimates = results.minimizer
CSV.write("../Output/fused_ridge_estimates_"*string(gamma)*".csv", DataFrame(estimates=estimates))

# Compute and export wage and rent beliefs
α = 0.24 
ϵ = estimates[1]
fe_i_ridge = estimates[2:num_i + 1]
fe_j_ridge = estimates[num_i+2:end]
rentbelief = coalesce.(exp.(-fe_i_ridge ./ (ϵ * α)), Inf)
wagebelief = coalesce.(exp.(fe_j_ridge ./ ϵ), 0.0)
CSV.write("../Output/fused_ridge_rentbelief_"*string(gamma)*".csv", DataFrame(wagebelief=rentbelief))
CSV.write("../Output/fused_ridge_wagebelief_"*string(gamma)*".csv", DataFrame(wagebelief=wagebelief))
