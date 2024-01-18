module SVD
using LinearAlgebra
using Statistics 
using Plots
using LaTeXStrings

include("../fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator, corr_noise_gen

include("../utils/paths.jl")

include("../fourier/FourierAnalysis.jl")
using .FourierAnalysis: intercept_and_exponent_from_log_psd

export create_ts_matrix, centralize_matrix, plot_eigen_spectrum

function filter_eigen_vals_array(atol::Float64,M::Matrix{Float64})
    prop_vals = svd(M).S  
    return filter(u -> u > atol, prop_vals)
end

"""
    create_ts_matrix(beta::Float64, number_of_observations=10::Int)::Matrix{Float64}

Returns Float64 a matrix by stacking a given number of observations and a linear correlation exponent
"""
function create_ts_matrix(beta::Float64, number_of_observations=10::Int)::Matrix{Float64}
    corr_noises_arr  = Float64[]
    corr_noises_arr = Array{Float64,1}[]

    for i in 1:number_of_observations
        corr_noise = CorrelatedNoise.correlated_noise_generator(1000,1.0,beta) #times series 
        push!(corr_noises_arr,corr_noise)
    end
    
    return M = hcat(corr_noises_arr...)'
end

function centralize_matrix(M::Matrix{Float64})::Matrix{Float64}
    mean_by_row = mean.(eachrow(M)) 
    return M .- mean_by_row
end

"""
   plot_eigen_spectrum(dir_to_save::String, M::Matrix{Float64})

Persist a pdf file with the plot_eigen_spectrum of the matrix M     
"""
function plot_eigen_spectrum(dir_to_save::String, M::Matrix{Float64})
    #build x, y axis 
    filtered_eig_vals = filter_eigen_vals_array(0.001,M)
    x = collect(Float64,1:length(filtered_eig_vals))
    
    #compute linear fit 
    params = FourierAnalysis.intercept_and_exponent_from_log_psd(x,filtered_eig_vals)
    
    #persist graph if doesn't exist
    full_file_path = joinpath(dir_to_save,"eigenspectrum.pdf")
    if !isfile(full_file_path)
        #plot styling
        plt = plot(x,filtered_eig_vals, label=L"{Eig val}_n", legend=false, xscale=:log10, yscale=:log10,alpha=0.2)
        #linear fit
        plot!(u -> exp10(params[1] + params[2]*log10(u)),minimum(x),maximum(x), xscale=:log10,yscale=:log10,lc=:red)
        
        title!("Eigen spectrum, beta_fit = $(round(-params[2],digits=3))")
        xlabel!(L"n")
        ylabel!("Eigen spectrum")
        
        #file saving
        savefig(plt, full_file_path)
    end
end
end #end of module





