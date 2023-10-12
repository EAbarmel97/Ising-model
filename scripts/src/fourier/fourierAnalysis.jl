module fourierAnalysis
export  compute_rfft, compute_psd, write_rfft,plot_psd

using FFTW
using Plots
using LaTeXStrings

include("../ising.jl")
using .ising: CRITICAL_TEMP

include("../utilities.jl")
using .utilities: get_array_from_txt

#= Auxiliray constants =#
const PSD_GRAPHS =  "graphs/psd"
const AUTOMATED_PSD_GRAPHS =  "graphs/automated/psd"

#= Function wrapper to the FFTW's rfft method. This wrapper takes a .txt file of numbers, parses them 
   an computes the rfft =#
function compute_rfft(file_path :: AbstractString) :: Array{ComplexF64,1}
    time_series = utilities.get_array_from_txt(file_path)
    return rfft(time_series) 
end 

#= Function to write over a .txt file a vector with the rfft of a signal(time series) =#
function write_rfft(arr :: Array{ComplexF64,1}, destination_dir :: AbstractString, 
    at_temp :: Float64, run :: Int)

    if at_temp != ising.CRITICAL_TEMP    
        rounded_temp = round(at_temp, digits=2)
        str_rounded_temp = replace("$rounded_temp","." => "_")
        file_name = "$destination_dir/rfft_global_magnetization_$(str_rounded_temp)_r$run.txt"
    else
        str_Tc_temp = replace("$at_temp","$(ising.CRITICAL_TEMP)" => "Tc",)
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
function compute_psd(arr :: Array{T,1} where T <: Complex) :: Array{Float64,1} 
    return abs2.(arr)
end

#= Function to determine the array containing sampling frecuencies when an abs path is given as argument =#
function sampling_freq_arr(file_path :: AbstractString) :: Array{Float64,1}
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

#= method to take the average psd when array of psd at different runs is given =#
function mean_psd(psd_array :: Array{Array{Float64,1},1}) :: Array{Float64,1}
    sum = zeros(length(psd_array[1]))
    for i in eachindex(psd_array)
        psd = psd_array[i]
        sum += psd    
    end
    return sum/length(psd_array)
end

function intercept_and_order_coef(x::Array{Float64,1},y::Array{Float64,1})::Array{Float64,1}
    X = hcat(ones(length(x)),x)
    return inv(X'*X)*(X'*y)
end
#= 
Module method for plotting psd wuth options to plot several psd on the same canvas, providing one generic 
under which all psd will be saved. 
    
NOTE: the paths to the rffts .txt files need to be absulute path/s 
=#
function plot_psd(temp_name_dir :: AbstractString, destination_dir :: AbstractString)
    
    #auxiliar variables
    curr_dir = pwd() 
    full_file_path = ""
    simuls_dir = "" 
    
    #if plots will be saved into graphs/automated then they were generated by noninteractively
    if contains(destination_dir,"automated")
        simuls_dir = joinpath(curr_dir,"all_simulations/automated")
    else     
        simuls_dir = joinpath(curr_dir,"all_simulations") 
    end
    
    dashed_str_temp = replace(temp_name_dir, "simulations_" => "")
    at_temp = joinpath(curr_dir,destination_dir,dashed_str_temp) # subdir ../graphs/automated/psd/T_x_y_z or ../graphs/psd/T_x_y_z
    mkpath(at_temp) #sub dir graphs/psd/T_x_y_z
    
    rffts_at_temp = joinpath(simuls_dir,temp_name_dir,"fourier")
    rffts_file_names = readdir(rffts_at_temp) #names of the rffts .txt file saved under the dir ../fourier/
    
    NUM_RUNS = length(readdir(rffts_at_temp))
    psd_array = Array{Float64,1}[]#array containing psd of different runs
    for run in 1:NUM_RUNS
        rfft_file_name = rffts_file_names[run]
        rfft_path = joinpath(rffts_at_temp,rfft_file_name)#abs path to the strigified file  for the rfft 

        #fetching and appending psd to the array containg the power spectra densities
        rfft = utilities.get_array_from_txt(Complex{Float64},rfft_path) #rfft of the M_n with initial temperature x_y_z
        deleteat!(rfft,(1,length(rfft))) #discarting the DC associated entry and the last element array 
        rfft = convert.(ComplexF64,rfft) #casting array to ComplexF64
        
        psd = fourierAnalysis.compute_psd(rfft) #array with the psd associated with RFFT[M_n]
        push!(psd_array,psd)
    end
    
    average_psd = mean_psd(psd_array) #mean psd

    #string manipulations
    magn_dir_at_temp = joinpath(simuls_dir,temp_name_dir,"magnetization")
    #all magnetization time series files have the number of lines, so the first file is picked up
    magn_file_name = readdir(magn_dir_at_temp)[1]
    magn_ts_abs_path = joinpath(magn_dir_at_temp,magn_file_name)
    
    #sampling frecuencies
    f = sampling_freq_arr(magn_ts_abs_path)

    log10_f = log10.(f)
    log10_mean_psd = log10.(average_psd)
    beta0, beta1 = intercept_and_order_coef(log10_f,log10_mean_psd)
    params = [beta0,beta1]

    #plot styling
    plt = plot(f, psd_array, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,alpha=0.2) #plot reference 
    
    plot!(f, average_psd, label=L"PSD \ \left( f \right)", legend=false, xscale=:log10, yscale=:log10,lc=:red)

    plot!((x) -> exp10(params[1] + params[2]*log10(x)),minimum(f),maximum(f),legend=false, xscale=:log10,yscale=:log10,lc=:black)
    
    str_temp = replace(dashed_str_temp,"T_" => "","_" => ".")
    title!("PSD for ts with init temp $(str_temp)")
    xlabel!(L"f")
    ylabel!("power density spectra")
    
    #file saving
    full_file_path = joinpath(at_temp,"psd_$(dashed_str_temp)_r1_$(NUM_RUNS).pdf")
    savefig(plt, full_file_path)
end
end #end of module 