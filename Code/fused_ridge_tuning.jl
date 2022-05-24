using Base.Threads, CSV, DataFrames, JLD, Optim


@time data = sort(DataFrame(load("../Data/nyc2010_lodes_wzero_wdelta.dta")), [:j, :i])
num_i = length(unique(data.i))
num_j = length(unique(data.j))
i_df = DataFrame(id_i=collect(1:num_i), i=sort(unique(data.i)))
j_df = DataFrame(id_j=collect(1:num_j), j=sort(unique(data.j)))
i_df.fe_i_ppml = CSV.read("../Data/ppml_fe_i.csv", DataFrame).fe_i_ppml
j_df.fe_j_ppml = CSV.read("../Data/ppml_fe_j.csv", DataFrame).fe_j_ppml


function calculate_MSE(ridge_df::DataFrame; data=data, num_i=num_i, num_j=num_j, i_df=i_df, j_df=j_df)
    ϵ_ridge = ridge_df.estimates[1]
    i_df.fe_i_ridge = ridge_df.estimates[2:num_i + 1]
    j_df.fe_j_ridge = ridge_df.estimates[num_i+2:end]
    data = outerjoin(data, i_df, on=:i);
    data = outerjoin(data, j_df, on=:j);
    X_ij_pred = exp.(ϵ_ridge .* data.log_delta + data.fe_i_ridge + data.fe_i_ridge)
    return sum(skipmissing((data.X_ij .- X_ij_pred) .^ 2)) / length(data.X_ij)
end


file_vec = Vector{String}()
γ_vec = Vector{Float64}()
MSE_vec = Vector{Float64}()
ϵ_ridge_vec = Vector{Float64}()
for file in readdir("../Output/")
    if occursin("fused_ridge_estimates", file)
        push!(file_vec, file)
        push!(γ_vec, parse(Float64,  split(split(file, "_")[end], ".csv")[1]))
        df = CSV.read("../Output/"*file, DataFrame)
        ϵ_ridge = df.estimates[1]; push!(ϵ_ridge_vec, ϵ_ridge);
        MSE = calculate_MSE(df)
        @time push!(MSE_vec, MSE)
        println(γ_vec[end], " ", ϵ_ridge, " ", MSE)
    end
end
df_tuning = sort(DataFrame(gamma=γ_vec, MSE=MSE_vec, epsilon_est=ϵ_ridge_vec), [:gamma])

min_idx = findmin(df_tuning.MSE[2:end])[2]
gamma_optimal = df_tuning.gamma[min_idx]
file_optimal = "../Output/fused_ridge_estimates_"*string(gamma_optimal)*".csv"
df_optimal = CSV.read(file_optimal, DataFrame)
CSV.write("../Output/fused_ridge_optimal.csv", df_optimal)

ϵ_ridge_optimal = df_optimal.estimates[1]
i_df.fe_i_ridge_optimal = df_optimal.estimates[2:num_i + 1]
j_df.fe_j_ridge_optimal = df_optimal.estimates[num_i+2:end]
CSV.write("../Output/fused_ridge_fe_i_optimal.csv", i_df)
CSV.write("../Output/fused_ridge_fe_j_optimal.csv", j_df)



