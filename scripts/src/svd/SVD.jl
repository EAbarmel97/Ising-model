module SVD
using LinearAlgebra
using Statistics 
using Base.Threads

using Plots
using LaTeXStrings

include("../fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator

include("../utils/paths.jl")

include("../fourier/FourierAnalysis.jl")
using .FourierAnalysis: intercept_and_exponent_from_log_psd

export create_ts_matrix, centralize_matrix, plot_eigen_spectrum, write_beta_beta_fit, read_beta_beta_fit, plot_beta_beta_fit, average_eigen_spectrum, BETA_BETA_FIT_FILE_PATH

const BETA_BETA_FIT_FILE_PATH  = joinpath(REPO_DIR, "beta_beta_fit.txt")

function pad_vector(array::Vector{T}, co_dim::Int64)::Vector{T} where {T <: Real}
    return vcat(array,zeros(T,co_dim))
end

function _fill_matrix_with_eigspectrum(eigen_spectra::Vector{Vector{T}}, max_rank::Int64)::Matrix{T} where {T <: Real}
    eigen_spectra_matrix = zeros(T, length(eigen_spectra), max_rank)
    
    for i in eachindex(eigen_spectra)
        eigen_spectra_matrix[i,:] .= pad_vector(eigen_spectra[i], max_rank - length(eigen_spectra[i]))
    end
    return eigen_spectra_matrix
end

"""
   average_eigen_spectrum(beta::Float64; number_of_realizations=100::Int,number_of_observations=1000::Int, num_samples=100::Int64)::Vector{Float64}

Gives the mean eigenspectrum averaged over a given number (num_samples) of eigen spectra of a matrix containing stacked realizations (number_of_realizations) 
of beta correlated time series with a fixed length (number_of_observations)
"""
function average_eigen_spectrum(beta::Float64; number_of_realizations=100::Int,number_of_observations=1000::Int,num_samples=100::Int64)
    eigen_spectra = Vector{Float64}[]
    ranks = Int64[]

    Threads.@threads for i in 1:num_samples
       eigen_spectrum = compute_eigenvals_from_beta_generated_noise(beta; number_of_realizations=number_of_realizations,number_of_observations=number_of_observations)
       push!(eigen_spectra, eigen_spectrum)
       push!(ranks, length(eigen_spectra[i]))
    end

    max_rank = maximum(ranks)
    eigen_spectra_matrix = _fill_matrix_with_eigspectrum(eigen_spectra, max_rank)
    cols = eachcol(eigen_spectra_matrix)
    
    return mean.(cols)
end

"""
    create_ts_matrix(beta::Float64, number_of_observations=10::Int)::Matrix{Float64}

Returns Float64 a matrix by stacking a given number of observations and a linear correlation exponent
"""
function create_ts_matrix(beta::Float64; number_of_realizations=100::Int,number_of_observations=1000::Int)::Matrix{Float64}
    corr_noises_arr  = Float64[]
    corr_noises_arr = Array{Float64,1}[]

    for _ in 1:number_of_realizations
        corr_noise = CorrelatedNoise.correlated_noise_generator(number_of_observations,1.0,beta) #times series 
        push!(corr_noises_arr,corr_noise)
    end
    
    return M = hcat(corr_noises_arr...)'
end

"""
   create_ploting_axes(M::Matrix{Float64})::Array{Array{Float64,1},1}

Returns an array of dim 2 whose elements are arrays of float64. Respectively they are the x-axis and y-exis to by ploted latter on
"""
function _create_ploting_axes(eigvals::Array{Float64})::Array{Array{Float64,1},1}
    axes = (collect(Float64,1:length(eigvals)),eigvals)
    return collect(axes)
end

function _create_ploting_axes(betas::Array{Float64}, fitted_betas::Array{Float64,1})::Array{Array{Float64,1},1}
    axes = (betas,fitted_betas)
    return collect(axes)
end

function _create_eigenspectrum_plot_file_path(dir_to_save,beta::Float64)
    file_name = replace("eigenspectrum_beta_$(round(beta,digits=3))", "." => "_") * ".pdf"
    return joinpath(dir_to_save,file_name)
end

"""
   create_beta_beta_fit_file()

Create a .txt file where the beta and beta fit will be written; When created it also writes the headers of the file
"""
function _create_beta_beta_fit_file()
    if !isfile(BETA_BETA_FIT_FILE_PATH)
        touch(BETA_BETA_FIT_FILE_PATH) 
        #writes headers to file one it is created
        open(BETA_BETA_FIT_FILE_PATH,"a+") do io
            write(io,"beta,beta_fit\n")
        end    
    end
end


function _write_beta_beta_fit_to_file(full_file_path::String,beta::Float64,beta_fit::Float64;is_eof::Bool=true)
    value_to_write = "$(beta),$(beta_fit)"
    if is_eof
        open(full_file_path,"a+") do io
            write(io, value_to_write)
        end 
    else
        open(full_file_path,"a+") do io
            write(io,value_to_write * "\n")
        end    
    end      
end


"""
    write_beta_beta_fit(from_beta::Float64, to_beta::Float64, num_of_betas::Int)

Writes over a .txt file the values of the several beta vs beta fit values. Th graphs of each individual beta fit can be persisted if wanted using the argument 
save_individual_plot. It's default value is false.
"""
function write_beta_beta_fit(from_beta::Float64, to_beta::Float64, num_of_betas::Int; number_of_realizations=10::Int, number_of_observations=1000::Int, num_samples=10::Int64)
    _create_beta_beta_fit_file()
    #= TO OPTIMIZE USING multi-threads =#
    Threads.@threads for beta in LinRange(from_beta, to_beta, num_of_betas)  
        is_eof = beta == to_beta #line will no be written to file if eof
        
        average_eigvals = average_eigen_spectrum(beta; 
                            number_of_realizations=number_of_realizations,
                            number_of_observations=number_of_observations,
                            num_samples=num_samples)                    

        params = compute_linear_fit_params(average_eigvals)

        _write_beta_beta_fit_to_file(BETA_BETA_FIT_FILE_PATH,beta,-params[2];is_eof=is_eof)
    end
end

"""
   read_beta_beta_fit(file_path::String)::Array{Array{Float64,1},1}

Reads teh contents of the .txt file with the values of the beta, beta_fit and returns them in a nested array
"""
function read_beta_beta_fit(file_path::String)::Array{Array{Float64,1},1}
    beta_array = Float64[] 
    beta_fit_array = Float64[]
    betas = (beta_array, beta_fit_array)

    lines = readlines(file_path)
    for line in lines[2:end]
        str_beta_beta_fit_arr = string.(split(line,","))
        beta_beta_fit = Base.parse.(Float64,str_beta_beta_fit_arr)
        push!(beta_array,beta_beta_fit[1])
        push!(beta_fit_array,beta_beta_fit[2])
    end

    return collect(betas)
end

"""
   compute_linear_fit_params(eigvals::Array{Float64,1})::Arraty{Float64,1}

Returns the parameters of the linear fit of an eigenspectrum 
"""
function compute_linear_fit_params(eigvals::Array{Float64,1})::Array{Float64,1}
    return FourierAnalysis.intercept_and_exponent_from_log_psd(collect(Float64,1:length(eigvals)),eigvals)
end

function centralize_matrix(M::Matrix{Float64})::Matrix{Float64}
    mean_by_row = mean.(eachrow(M)) 
    return M .- mean_by_row
end

"""
   filter_singular_vals_array(atol::Float64,M::Matrix{Float64})

Returns an array of singular valuesfiltered by absulte tolerance
"""
function filter_singular_vals_array(M::Matrix{Float64};atol=eps(Float64)::Float64)
    prop_vals = svd(M).S  
    return filter(u -> u > atol, prop_vals)
end

function compute_eigvals(M::Matrix{Float64}; drop_first=true::Bool)::Array{Float64,1}
    if drop_first 
        return abs2.(filter_singular_vals_array(M))[2:end]
    end
    
    return abs2.(filter_singular_vals_array(M))
end

function compute_eigenvals_from_beta_generated_noise(beta::Float64; number_of_realizations=10::Int,number_of_observations=1000::Int)::Array{Float64,1}
    ts_matrix = create_ts_matrix(beta;number_of_realizations=number_of_realizations,number_of_observations=number_of_observations)
    M = centralize_matrix(ts_matrix)
    eigvals = compute_eigvals(M)

    return eigvals
end

function compute_eigenspectrum_length(beta::Float64; number_of_realizations=10::Int,number_of_observations=1000::Int)::Int64
    eigvals = compute_eigenvals_from_beta_generated_noise(beta; number_of_realizations=number_of_realizations,number_of_observations=number_of_observations)
    return length(eigvals)
end

"""
   plot_eigen_spectrum(dir_to_save::String, M::Matrix{Float64})

Persist a pdf file with the plot_eigen_spectrum of the matrix M     
"""
function _plot_eigen_spectrum(dir_to_save::String, eigvals::Array{Float64,1},beta::Float64)
    #build x, y axis; y being the eigenspectrum and x it's enumeration
    ploting_axes = _create_ploting_axes(eigvals)

    #compute linear fit 
    params = compute_linear_fit_params(ploting_axes[2])
    
    full_file_path = _create_eigenspectrum_plot_file_path(dir_to_save,beta)
    
    #persist graph if doesn't exist
    if !isfile(full_file_path)
        #plot styling
        plt = plot(ploting_axes[1],ploting_axes[2], label=L"{Eig val}_n", legend=false, xscale=:log10, yscale=:log10,alpha=0.2)
        #linear fit
        plot!(u -> exp10(params[1] + params[2]*log10(u)),minimum(ploting_axes[1]),maximum(ploting_axes[1]), xscale=:log10,yscale=:log10,lc=:red)
        
        title!("Eigen spectrum,beta = $(round(beta,digits=3)), beta_fit = $(round(-params[2],digits=3))")
        xlabel!(L"n")
        ylabel!("Eigen spectrum")
        
        #file saving
        savefig(plt, full_file_path)
    end
end

function plot_eigen_spectrum(dir_to_save::String,M::Matrix{Float64},beta::Float64)
    eigvals = compute_eigvals(M)
    _plot_eigen_spectrum(dir_to_save,eigvals,beta)
end


"""
   plot_beta_beta_fit(file_path::String)

Persists a graph of the beta vs beta fit is a given dir
"""
function plot_beta_beta_fit(file_path::String)

    beta_array, beta_fit_array = read_beta_beta_fit(file_path)
    
    beta_beta_fit_plot_file_path = joinpath(AUTOMATED_EIGEN_SEPCTRUM_GRAPHS_DIR,"beta_vs_beta_fit.pdf")

    if !isfile(beta_beta_fit_plot_file_path)
        #plot styling
        x = collect(1:length(beta_array))
        plt = plot(x,[beta_array, beta_fit_array],label=["beta" "beta_fit"], seriestype=:scatter,alpha=0.5)
        #linear fit
        plot!(u -> u, minimum(x), maximum(x),label="id",lc=:black)
        
        title!("beta vs beta_fit, beta from $(beta_array[1]) to $(beta_array[end])")

        #file saving
        savefig(plt, beta_beta_fit_plot_file_path)
    end
end

"""
   plot_beta_beta_fit(file_path::String)

Persists a graph of the beta vs beta fit is a given dir
"""
function plot_beta_beta_fit(func::Function; file_path::String=BETA_BETA_FIT_FILE_PATH)
    axes = _create_ploting_axes(func(file_path)...)

    beta_beta_fit_plot_file_path = joinpath(AUTOMATED_EIGEN_SEPCTRUM_GRAPHS_DIR,"beta_vs_beta_fit.pdf")

    if !isfile(beta_beta_fit_plot_file_path)
        #plot styling
        plt = plot(axes[1],axes[2],label=["beta" "beta_fit"], seriestype=:scatter,alpha=0.5)
        #linear fit
        plot!(u -> u,label="id",lc=:black)
        
        title!("beta vs beta_fit, beta from $(axes[1][1]) to $(axes[1][end])")

        #file saving
        savefig(plt, beta_beta_fit_plot_file_path)
    end
end

end #end of module





