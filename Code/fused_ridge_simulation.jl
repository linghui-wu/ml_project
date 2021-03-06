using Base.Threads, CSV, DataFrames, JLD, Optim

dgp = ARGS[1]
event = parse(Int64, ARGS[2])

data = sort(CSV.read("../Data/simulated_data/"*dgp*"_"*string(event)*".csv", DataFrame), [:n, :k])
num_k = length(unique(data.k))
num_n = length(unique(data.n))
k_df = DataFrame(id_k=collect(1:num_k), k=sort(unique(data.k)))
n_df = DataFrame(id_n=collect(1:num_n), n=sort(unique(data.n)))

fe_k = unique(DataFrame(n=data.n, origin_fe=data.origin_fe)).origin_fe
fe_n = unique(DataFrame(k=data.k, dest_fe=data.dest_fe)).dest_fe

neighboring_tracts = CSV.read("../Data/simulated_data/tract_adjacency.csv", DataFrame)
rename!(neighboring_tracts, :source=>:k)
k_neighboring = innerjoin(k_df, neighboring_tracts, on=:k)

neighboring_tracts = CSV.read("../Data/simulated_data/tract_adjacency.csv", DataFrame)
rename!(neighboring_tracts, :source=>:n)
n_neighboring = innerjoin(n_df, neighboring_tracts, on=:n)

ϵ = -8.0
params = vcat(ϵ, fe_k, fe_n)

function ll(params; γ=gamma, data=data, num_k=num_k, num_n=num_n, D_k=k_neighboring, D_n=n_neighboring)
    ϵ = params[1]
    α_k_vec = params[2:num_k + 1]
    α_n_vec = params[num_k+2:end]
    k_df = DataFrame(id_k=collect(1:num_k), k=unique(data.k))
    k_df.fe_k = α_k_vec 
    n_df = DataFrame(id_n=collect(1:num_n), n=unique(data.n))
    n_df.fe_n = α_n_vec 
    data = outerjoin(data, k_df, on=:k)
    data = outerjoin(data, n_df, on=:n)
    term1 = -data.ell_kn .* (ϵ .* log.(data.delta) .+ data.fe_k .+ data.fe_n)
    term2 = data.delta .^ ϵ .* exp.(data.fe_k .+ data.fe_n)
    k_nbhd = outerjoin(k_df, D_k, on=:k, makeunique=true)
    rename!(k_nbhd, :k=>:k_origin, :fe_k=>:fe_k_origin, :neighbor=>:k)
    k_nbhd = outerjoin(k_df, k_nbhd, on=:k, makeunique=true, matchmissing=:equal)
    k_nbhd = k_nbhd[completecases(k_nbhd), :]  # Drop row with missing values in columns
    n_nbhd = outerjoin(n_df, D_n, on=:n, makeunique=true)
    rename!(n_nbhd, :n=>:n_origin, :fe_n=>:fe_n_origin, :neighbor=>:n)
    n_nbhd = outerjoin(n_df, n_nbhd, on=:n, makeunique=true, matchmissing=:equal)
    n_nbhd = n_nbhd[completecases(n_nbhd), :]  # Drop row with missing values in columns
    penalty_k = sum((k_nbhd.fe_k .- k_nbhd.fe_k_origin) .^ 2)
    penalty_n = sum((n_nbhd.fe_n .- n_nbhd.fe_n_origin) .^ 2)
    log_likelihood_val = sum(skipmissing(term1) .+ skipmissing(term2)) + γ * (penalty_k + penalty_n)
    # println(ϵ, " ", log_likelihood_val)
    return log_likelihood_val
end

function g!(G, params; data=data, γ=gamma, num_k=num_k, num_n=num_n, D_k=k_neighboring, D_n=n_neighboring)
    ϵ = params[1]
    α_k_vec = params[2:num_k + 1]
    α_n_vec = params[num_k+2:end]
    k_df = DataFrame(id_k=collect(1:num_k), k=unique(data.k))
    k_df.fe_k = α_k_vec 
    n_df = DataFrame(id_n=collect(1:num_n), n=unique(data.n))
    n_df.fe_n = α_n_vec 
    data = outerjoin(data, k_df, on=:k)
    data = outerjoin(data, n_df, on=:n)
    k_nbhd = outerjoin(k_df, D_k, on=:k, makeunique=true)
    rename!(k_nbhd, :k=>:k_origin, :fe_k=>:fe_k_origin, :neighbor=>:k)
    k_nbhd = outerjoin(k_df, k_nbhd, on=:k, makeunique=true, matchmissing=:equal)
    k_nbhd.diff_term = 2 * 2 * γ * (k_nbhd.fe_k_origin .- k_nbhd.fe_k)
    k_nbhd = k_nbhd[completecases(k_nbhd), :]  # Drop row with missing values in columns
    n_nbhd = outerjoin(n_df, D_n, on=:n, makeunique=true)
    rename!(n_nbhd, :n=>:n_origin, :fe_n=>:fe_n_origin, :neighbor=>:n)
    n_nbhd = outerjoin(n_df, n_nbhd, on=:n, makeunique=true, matchmissing=:equal)
    n_nbhd.diff_term = 2 * γ * (n_nbhd.fe_n_origin .- n_nbhd.fe_n)
    n_nbhd = n_nbhd[completecases(n_nbhd), :] 
    G[1] = sum(skipmissing(log.(data.delta) .* (data.ell_kn .+ exp.(data.fe_k .+ data.fe_n) .* data.delta .^ ϵ)))
    Threads.@threads for idx in 1:num_k
        data_k = data[data.k .== idx, :];
        term1 = sum(skipmissing(data_k.ell_kn + exp.(data_k.fe_k .+ data_k.fe_n) .* data_k.delta .^ ϵ))
        k_nbhd_k = k_nbhd[k_nbhd.k_origin .== idx, :]
        term2 = sum(k_nbhd_k.diff_term)
        G[1 + idx] = (term1 + term2)
    end
    Threads.@threads for idx in 1:num_n
        data_n = data[data.n .== idx, :];
        term1 = sum(skipmissing(data_n.ell_kn + exp.(data_n.fe_k .+ data_n.fe_n).* data_n.delta .^ ϵ))
        n_nbhd_n = n_nbhd[n_nbhd.n_origin .== idx, :]
        term2 = sum(n_nbhd_n.diff_term)
        G[1 + num_k + idx] = (term1 + term2)        
    end
    return G
end


function calculate_MSE(ridge_df::DataFrame; data=data, num_k=num_k, num_n=num_n, k_df=k_df, n_df=n_df)
    ϵ_ridge = ridge_df.estimates[1]
    k_df.fe_k_ridge = ridge_df.estimates[2:num_k + 1]
    n_df.fe_n_ridge = ridge_df.estimates[num_n+2:end]
    data = outerjoin(data, k_df, on=:k);
    data = outerjoin(data, n_df, on=:n);
    ell_kn_pred = exp.(-ϵ_ridge .* log.(data.delta) .+ data.fe_n_ridge .+ data.fe_k_ridge)
    return sum(skipmissing((data.ell_kn .- ell_kn_pred) .^ 2)) / length(data.ell_kn)
end

gamma_vec = collect(0:50:1500)
MSE_vec = zeros(length(gamma_vec))
for idx in collect(1:length(gamma_vec))
    gamma = gamma_vec[idx]

    lower = -15 * ones(num_k + num_n + 1); lower[1] = -8.5;
    upper = 15 * ones(num_k + num_n + 1); upper[1] = -7.5;
    @time results = optimize(ll, g!, lower, upper, params, Fminbox(GradientDescent()), 
            Optim.Options(show_trace=true, x_abstol=1e-1, f_tol=1.0, g_tol=1.0, 
                iterations=2, outer_iterations=2))
    estimates_gamma = results.minimizer
    CSV.write("../Output_simulation/fused_ridge_simulation_estimates_"*dgp*"_nbhd_"*string(event)*"_"*string(gamma)*".csv", 
        DataFrame(estimates=estimates_gamma))
    MSE_vec[idx] = calculate_MSE(DataFrame(estimates=estimates_gamma))
end
