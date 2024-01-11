
using LinearAlgebra
using Statistics 
using Plots

include("../src/fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator, corr_noise_gen

include("../src/utils/paths.jl")

corr_noises_arr  = Float64[]
corr_noises_arr = Array{Float64,1}[]

for i in 1:10
    beta_rn = rand()
    corr_noise = CorrelatedNoise.correlated_noise_generator(1000,1.0,beta_rn) #times series 
    push!(corr_noises_arr,corr_noise)
end    

M = hcat(corr_noises_arr...)' #matrix of time series produced from correlated noises 
@show size(M)
@show row_dim,col_dim = size(M)[1],size(M)[2]

mean_by_observation = reshape(mean.(eachcol(M)),col_dim,1) #row vector with per-row mean 
ones = reshape(fill(1,row_dim),row_dim,1) #col vector with ones 

@show size(mean_by_observation), size(ones)
centered_M = M - ones * mean_by_observation'

prop_vals = svd(centered_M).S

len = length(prop_vals)
plt = plot!(x=collect(1:len),prop_vals)
savefig(plt, "$(pwd())/scripts/svd/eigvals.pdf")




