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

function calculate_MSE(ridge_df::DataFrame; data=data, num_k=num_k, num_n=num_n, k_df=k_df, n_df=n_df)
    ϵ_ridge = ridge_df.estimates[1]
    k_df.fe_k_ridge = ridge_df.estimates[2:num_k + 1]
    n_df.fe_n_ridge = ridge_df.estimates[num_n+2:end]
    data = outerjoin(data, k_df, on=:k);
    data = outerjoin(data, n_df, on=:n);
    ell_kn_pred = exp.(ϵ_ridge .* log.(data.delta) .+ data.fe_n_ridge .+ data.fe_k_ridge)
    return sum(skipmissing((data.ell_kn .- ell_kn_pred) .^ 2)) / length(data.ell_kn)
end

file_vec = Vector{String}()
γ_vec = Vector{Float64}()
MSE_vec = Vector{Float64}()
ϵ_ridge_vec = Vector{Float64}()
for file in readdir("../Output/")
    if occursin("fused_ridge_simulation_estimates_"*dgp*"_nbhd_"*string(event)*"_", file)
        push!(file_vec, file)
        push!(γ_vec, parse(Float64,  split(split(file, "_")[end], ".csv")[1]))
        df = CSV.read("../Output/"*file, DataFrame)
        ϵ_ridge = df.estimates[1]; push!(ϵ_ridge_vec, ϵ_ridge);
        push!(MSE_vec, calculate_MSE(df))
    end
end
df_tuning = sort(DataFrame(gamma=γ_vec, MSE=MSE_vec, epsilon_est=ϵ_ridge_vec), [:gamma])
min_idx = findmin(df_tuning.MSE)[2]
gamma_optimal = df_tuning.gamma[min_idx]
file_optimal = "../Output/fused_ridge_simulation_estimates_"*dgp*"_nbhd_"*string(event)*"_"*string(gamma_optimal)*".csv"
df_optimal = CSV.read(file_optimal, DataFrame)
CSV.write("../Output/fused_ridge_simulation_"*dgp*"_optimal.csv", df_optimal)

ϵ_ridge_optimal = df_optimal.estimates[1]
k_df.fe_k_ridge_optimal = df_optimal.estimates[2:num_k + 1]
n_df.fe_n_ridge_optimal = df_optimal.estimates[num_n+2:end]
CSV.write("../Output/fused_ridge_simulation_"*dgp*"_fe_i_optimal.csv", k_df)
CSV.write("../Output/fused_ridge_simulation_"*dgp*"_fe_j_optimal.csv", n_df)
