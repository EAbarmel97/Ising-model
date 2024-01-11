
using LinearAlgebra
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

@show M = hcat(corr_noises_arr...)' #matrix of time series produced from correlated noises 
@show rank(M)

prop_vals = svd(M).S
len = length(prop_vals)

plt = plot!(x=collect(1:len),prop_vals)
savefig(plt, "$(pwd())/scripts/svd/eigvals.pdf")





