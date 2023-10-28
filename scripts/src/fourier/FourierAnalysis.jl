module FourierAnalysis
export  compute_rfft, compute_psd, write_rfft, plot_psds

using FFTW
using Plots
using LaTeXStrings

include("../Ising.jl")
using .Ising: CRITICAL_TEMP

include("../utils/utilities.jl")
using .utilities: get_array_from_txt

include("../utils/paths.jl")

"""
    compute_rfft(file_path::AbstractString)::Array{ComplexF64,1}

Function wrapper to the FFTW's rfft method. Takes a .txt file of numbers, parses them 
an computes its rfft
"""
function compute_rfft(file_path::AbstractString)::Array{ComplexF64,1}
    time_series = utilities.get_array_from_txt(file_path)
    return rfft(time_series) 
end 

#= Function to write over a .txt file a vector with the rfft of a signal(time series) =#
function write_rfft(arr::Array{ComplexF64,1}, destination_dir::String, at_temp::Float64, run::Int64)

    if at_temp != Ising.CRITICAL_TEMP    
        rounded_temp = round(at_temp, digits=2)
        str_rounded_temp = replace("$rounded_temp","." => "_")
        file_name = "$destination_dir/rfft_global_magnetization_$(str_rounded_temp)_r$run.txt"
    else
        str_Tc_temp = replace("$at_temp","$(Ising.CRITICAL_TEMP)" => "Tc",)
        file_name = "$destination_dir/rfft_global_magnetization_$(str_Tc_temp)_r$run.txt"
    end

    #if file doesn't exist 
    if !isfile(file_name)    
        touch(file_name)
        fourier_transform_file = open(file_name,"a+")
        for i in eachindex(arr)
            write(fourier_transform_file,"$(arr[i])\n")
            #for the last index
            if i == length(arr)
               write(fourier_transform_file,"$(arr[i])")
            end    
        end
        close(fourier_transform_file)    
    end
end

#= Function to compute the power spectral density =#
function compute_psd(arr::Array{T,1}) where T <: Complex
    squared_norms = Array{Float64}[]
    squared_norms = abs2.(arr)
    return squared_norms
end

"""
    sampling_freq_arr(file_path::String)::Array{Float64,1}

Determines the array containing sampling frecuencies when an abs path is given as argument
"""
function sampling_freq_arr(file_path::String)::Array{Float64,1}
    freq_arr = []
    open(file_path,"r") do io
        num_lines = countlines(io)
        #converting to Array{Float64} to be able to use deleteat!
        sampling_freq_arr = rfftfreq(num_lines)
        freq_arr = convert.(Float64,sampling_freq_arr) 

        deleteat!(freq_arr,1)
    end    
    return freq_arr
end

function sampling_freq_arr(N::Int64)
    freq_arr = rfftfreq(N)
    freq_arr = convert.(Float64,freq_arr) 

    deleteat!(freq_arr,1)

    return freq_arr
end

""""
    mean_psd(psd_array::Array{Array{Float64,1},1})::Array{Float64,1}

returns the average psd when array of psd at different runs is given
"""
function mean_psd(psd_array::Array{Array{Float64,1},1})::Array{Float64,1}
    sum = zeros(length(first(psd_array)))
    for i in eachindex(psd_array)
        psd = psd_array[i]
        sum += psd    
    end

    return sum/length(psd_array)
end

#= Custom IO auxiliar functions =#

"""
auxiliar function to determine if the psd plot will be saved under /graphs/automated/psd or /graphs/psd/
"""
function determine_simulation_dir(destination_dir::String)::String
    if contains(destination_dir,"automated")
        simuls_dir = joinpath("all_simulations","automated")
    else     
        simuls_dir = "all_simulations"
    end

    return simuls_dir
end

"""
    create_graphs_temp_sub_dir(temp_name_dir::String,destination_dir::String)::Nothing

creates under /graphs/ the folder corresponding to a the simulations made a given temperature

Ex: /graphs/psd/T_1_34
"""
function create_graphs_temp_sub_dir(temp_name_dir::String,destination_dir::String)::Nothing
    dashed_str_temp = replace(temp_name_dir, "simulations_" => "")
    at_temp = joinpath(destination_dir,dashed_str_temp) # subdir ../graphs/automated/psd/T_x_y_z or ../graphs/psd/T_x_y_z
    mkpath(at_temp) #sub dir graphs/psd/T_x_y_z

    return nothing
end

"""
    file_names_in_fourier_dir(temp_name_dir::String, simuls_dir::String)::Array{String,1}

returns a list with all the rffts saved under a simulation dir at a given temperature
"""
function file_names_in_fourier_dir(temp_name_dir::String, simuls_dir::String)::Array{String,1}
    rffts_at_temp = joinpath(simuls_dir,temp_name_dir,"fourier")
    return readdir(rffts_at_temp) #names of the rffts .txt file saved under the dir ../fourier/
end

"""
    num_runs(temp_dir_name::String,simuls_dir)::Int64

computes the number of files contained in the ../fourier/ dir. That is the same as the number of runs
"""
function num_runs(temp_dir_name::String,simuls_dir)::Int64
    return length(file_names_in_fourier_dir(temp_dir_name,simuls_dir))
end

"""
    psd_arr_by_run(temp_dir_name::String,simuls_dir::String)::Array{Array{Float64,1},1}

returns an array with the psd at each run.
"""
function psd_arr_by_run(temp_dir_name::String,simuls_dir::String)::Array{Array{Float64,1},1}
    rffts_at_temp = joinpath(simuls_dir,temp_dir_name,"fourier")
    psd_array = Array{Float64,1}[]
    for rfft_file_name in file_names_in_fourier_dir(temp_dir_name,simuls_dir)
        rfft_path = joinpath(rffts_at_temp,rfft_file_name)#abs path to the strigified file  for the rfft 

        #fetching and appending psd to the array containg the power spectra densities
        rfft = utilities.get_array_from_txt(Complex{Float64},rfft_path)
        #rfft of the M_n with initial temperature x_y_z
        rfft = rfft[2:end-1]#discarting the DC associated entry and the last element array 
        
        rfft = convert.(ComplexF64,rfft) #casting array to ComplexF64
        
        psd = compute_psd(rfft) #array with the psd associated with RFFT[M_n]
        push!(psd_array,psd)
    end
    
    return psd_array
end

"""
    psd_graph_file_path(temp_dir_name::String,destination_dir::String)::String

Builds the absolute path to the psd graph. Of all superimposed psd ar different runs

Ex: "graphs/automated/psd_T_0_32_r_1_10.pdf"
"""
function psd_graph_file_path(temp_dir_name::String,destination_dir::String)::String
    NUM_RUNS = num_runs(temp_dir_name, determine_simulation_dir(destination_dir))
    
    t_dashed_str_temp = replace(temp_dir_name,"simulations_" => "")
    at_temp = joinpath(destination_dir,t_dashed_str_temp)
    full_file_path = joinpath(at_temp,"psd_$(t_dashed_str_temp)_r1_$(NUM_RUNS).pdf")
    
    return full_file_path
end

#= Auxiliar functions related to linear fit =#
"""
    determines_noise_or_movement(beta1::Float64)::String

returns a string corresponding to the classification of the simulated correlated noise or movement depending on the 
values of the beta1 linear fit parameter 
"""
function determines_noise_or_movement(beta1::Float64)::String
    if beta1 == 0.0
        description = "white noise"
    end
    
    if 0.0 < beta1 < 1.0
        description = "fractal noise with beta = $beta1"
    end     
    
    if beta1 == 1.0
        description = "pink noise"
    end

    if 1.0 < beta1 < 2.0
        description = "fractal brownian motion with beta = $beta1"
    end   

    if beta1 == 2.0
        description = "brownian motion"
    end
    
    return description
end

"""
    psd_graph_file_path(destination_dir::String,beta0::Float64,beta1::Float64)::String

Outputs the file path where the psd associated of the simulated time series will be saved 
"""
function psd_graph_file_path(destination_dir::String,A::Float64,beta::Float64)::String

    function separate_by_dashes(str::String)::String
        return replace(str, "." => "_")
    end
    
    full_file_path = joinpath(destination_dir,"psd_log_A_$(separate_by_dashes(string(round(log10(A),digits=2))))_beta_$(separate_by_dashes(string(round(beta,digits=2)))).pdf")

    return full_file_path
end


"""
    intercept_and_exponent(x::Array{Float64,1},y::Array{Float64,1})::Array{Float64,1}

Returns an array with the parameter estimators for a linear fit
"""
function intercept_and_exponent(x::Array{Float64,1},y::Array{Float64,1})::Array{Float64,1}
    X = hcat(ones(length(x)),x)

    return inv(X'*X)*(X'*y)
end

"""
    intercept_and_exponent_from_log_psd(f::Array{Float64,1},average_psd::Array{Float64,1})::Array{Float64,1}

Gives a 2 dimensional array containing the parameter estimators of a 2d linear fit
"""
function intercept_and_exponent_from_log_psd(f::Array{Float64,1},average_psd::Array{Float64,1})::Array{Float64,1}
    log10_f = log10.(f)
    log10_mean_psd = log10.(average_psd)
    beta0, beta1 = intercept_and_exponent(log10_f,log10_mean_psd)

    return [beta0,beta1]
end

function create_exponent_dir_and_file()
    exponent_dir = joinpath("graphs","automated","exponent")
    exponent_file_path = joinpath(exponent_dir,"exponent.txt")

    if !isdir(exponent_dir)
        mkpath(exponent_dir_dir)
        touch(exponent_file_path)
    end
end

function write_exponent(exponent::Float64,num_runs::Int64)
    exponent_file_path = "graphs/automated/exponent/exponent.txt"
    if i != num_runs
        val_to_append = string(exponent,"\n") 
    else
        val_to_append = string(exponent)
    end

    open(exponent_file_path,"a+") do exponent_file_path
        write(exponent_file_path,val_to_append)
    end       
end

"""

"""
function plot_log_PSD_exponent(log_PSD_exponent_file_path::String)

    if !isfile(log_PSD_exponent_file_path)

        f = sampling_freq_arr()
        #plot styling
        plt = plot(f, psd_array, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,alpha=0.2) #plot reference 
        
        plot!(f, average_psd, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,lc=:red)
        #linear fit
        plot!((x) -> exp10(params[1] + params[2]*log10(x)),minimum(f),maximum(f),legend=false, xscale=:log10,yscale=:log10,lc=:black)
        
    
        title!("PSD for ts with init temp $(str_temp)")
        xlabel!(L"f")
        ylabel!("power density spectra")
        
        #file saving
        savefig(plt, full_file_path)
    end
end

"""
    plot_psd(temp_name_dir::AbstractString,destination_dir::AbstractString)

Plots all psd in log-log superimposed on a same canvas, highlighting the mean psd in red, and the linear fit as well
"""
function plot_psd(temp_name_dir::AbstractString,destination_dir::AbstractString)
    simuls_dir = determine_simulation_dir(destination_dir)

    create_graphs_temp_sub_dir(temp_name_dir,destination_dir)

    psd_array = psd_arr_by_run(temp_name_dir,simuls_dir)
    average_psd = mean_psd(psd_array) #mean psd

    #string manipulations
    magn_dir_at_temp = joinpath(simuls_dir,temp_name_dir,"magnetization")
    #all magnetization time series files have the number of lines, so the first file is picked up
    magn_file_name = readdir(magn_dir_at_temp)[1]
    magn_ts_abs_path = joinpath(magn_dir_at_temp,magn_file_name)
    
    f = sampling_freq_arr(magn_ts_abs_path)
    params = intercept_and_exponent_from_log_psd(f,average_psd)

    full_file_path = psd_graph_file_path(temp_name_dir,destination_dir)

    if !isfile(full_file_path)
        #plot styling
        plt = plot(f, psd_array, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,alpha=0.2) #plot reference 
        #
        plot!(f, average_psd, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,lc=:red)
        #linear fit
        plot!((x) -> exp10(params[1] + params[2]*log10(x)),minimum(f),maximum(f),legend=false, xscale=:log10,yscale=:log10,lc=:black)
        
        str_temp = replace(temp_name_dir,"simulations_T_" => "", "_" => ".")
    
        title!("PSD for ts with init temp $(str_temp)")
        xlabel!(L"f")
        ylabel!("power density spectra")
        
        #file saving
        savefig(plt, full_file_path)
    end        
end

"""
    plot_psd(destination_dir::AbstractString)

Plots the psd associated with a complex array of numbers in log-log scale, highlighting the linear fit in black
"""
function plot_psd(x::Array{Float64,1}, destination_dir::String,A::Float64,beta::Float64)
    rfft_arr = rfft(x)
    f = sampling_freq_arr(length(x))

    psd = compute_psd(rfft_arr)
    deleteat!(psd,1)

   #=  @show any(iszero,psd)
    @show length(filter(iszero,psd))
    @show findall(x-> x==0,psd) =#

    params = intercept_and_exponent_from_log_psd(f,psd)
    
    full_file_path = psd_graph_file_path(destination_dir,A,beta)
    if !isfile(full_file_path)
        #plot styling
        plt = plot(f,psd, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,alpha=0.2) #plot reference 
        
        #expected linear fit
        #plot!((u) -> exp10(beta0-beta1*log10(u)),minimum(f),maximum(f),legend=false, xscale=:log10,yscale=:log10,lc=:black)
        plot!(f,(u) -> exp10(log10(A)-beta*log10(u)),minimum(f),maximum(f),legend=false, xscale=:log10,yscale=:log10,lc=:black)
        #linear fit
        plot!((u) -> exp10(params[1] + params[2]*log10(u)),minimum(f),maximum(f),legend=false, xscale=:log10,yscale=:log10,lc=:red)
        
        title!("PSD for ts with A = $A and beta = $beta")
        xlabel!(L"f")
        ylabel!("power density spectra")
        
        #file saving
        savefig(plt, full_file_path)
    end    
end
end #end of module 