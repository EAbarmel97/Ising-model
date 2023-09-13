module fourierAnalysis
export  compute_rfft, compute_psd, write_rfft,plot_psd

using FFTW
using Plots
using LaTeXStrings

include("../utilities.jl")
using .utilities: get_array_from_txt, replace_with_dict,neglect_N_first_from_array!
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
    rounded_temp = round(at_temp, digits=2)
    str_rounded_temp = replace("$rounded_temp","." => "_")
    file_name = "$destination_dir/rfft_global_magnetization_$(str_rounded_temp)_$run.txt"

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

#= Function to determine the array containing sampling frecuencies when an array of abs paths are given as argument =#
function sampling_freq_arr(file_paths :: AbstractArray) :: Array{Float64,1}
    aux_file_path = file_paths[1]
    freq_arr = []

    open(aux_file_path,"r+") do io
        num_lines = countlines(io)
        #converting to Array{Float64} to be able to use deleteat!
        sampling_freq_arr = rfftfreq(num_lines)
        freq_arr = convert.(Float64,sampling_freq_arr) 

        deleteat!(freq_arr,1)
    end     
    return freq_arr
end

#= function that returns the stringified temperature but with "_" as a separator =#
function get_str_dashed_temp(str_rfft_path :: AbstractString) :: AbstractString
    curr_dir = pwd()
    #replacement dictionary
    replace_dict = Dict(curr_dir => "","/all_simulations/automated/simulations_T_" => "",
    r"/fourier/rfft_global_magnetization_\d{1}_\d{1,2}_\d{1}.txt" => "")

    if contains(str_rfft_path,r"/all_simulations/automated/")
        str_dashed_temp = utilities.replace_with_dict(str_rfft_path, replace_dict)
    else
        str_dashed_temp = utilities.replace_with_dict(str_rfft_path, replace_dict)
    end

    return str_dashed_temp
end

#= 
Module method for plotting psd wuth options to plot a psd providing one generic 
under which all psd will be saved. 

NOTE: the path to the rfft .txt file need to be absulute 
=#

function plot_psd(str_rfft_path :: AbstractString, destination_dir :: AbstractString, run :: Int)
    #auxiliar defs
    curr_dir = pwd()
    full_file_path = destination_dir
    rel_path_sub_dir = ["magnetization","global_magnetization_r$run.txt"]

    #= abs path to the .txt file with the magnetization time series =# 
    sub_dirs_array = splitpath(str_rfft_path) #array containing all subdirs in str_rfft_path
    deleteat!(sub_dirs_array,length(sub_dirs_array) - 1 : length(sub_dirs_array))
    #array containing all subdirs included in the abs path the .txt with the magnetization time series
    append!(sub_dirs_array,rel_path_sub_dir) 
    magn_ts_abs_path = joinpath(sub_dirs_array) #path 

    #array of sampling frequencies
    f = sampling_freq_arr(magn_ts_abs_path)
    
    rfft = utilities.get_array_from_txt(Complex{Float64},str_rfft_path) #rfft of the M_n with initial temperature x_y_z
    deleteat!(rfft,(1,length(rfft))) #discarting the DC associated entry and the last element array 
   
    rfft = convert.(ComplexF64,rfft) #casting array to ComplexF64
    
    psd = fourierAnalysis.compute_psd(rfft) #array with the psd THE RFFT[M_n]
   
    #string manipulations
    str_dashed_temp = get_str_dashed_temp(str_rfft_path)
    full_file_path *= "psd_$(str_dashed_temp)_r$(run).pdf" 
    replace_dict = Dict("_" => ".")
    str_temp = replace_with_dict(str_dashed_temp,replace_dict)
    
    #= plot styling =#
    plt = plot(f, psd, label=L"PSD \ \left( f \right)") #plot reference 
    title!("PSD associated to a ts with init temp = $(str_temp)")
    xlabel!(L"f")
    ylabel!("power density spectrum")

    #= file saving  =#
    savefig(plt, full_file_path)
end

#= 
Module method for plotting psd wuth options to plot several psd on the same canvas, providing one generic 
under which all psd will be saved. 

NOTE: the paths to the rffts .txt files need to be absulute path/s 
=#
function plot_psd(str_rfft_paths :: AbstractArray, destination_dir :: AbstractString, run :: Int )
    curr_dir = pwd()
    full_file_path = destination_dir
    #auxiliar arrays
    str_dashed_temp_array = []
    psd_array = []

    
    for i in eachindex(str_rfft_paths)
        str_rfft_path = str_rfft_paths[i]

        #= fetching and appending psd to the array containg the power spectra densities =#
        rfft = utilities.get_array_from_txt(Complex{Float64},str_rfft_path) #rfft of the M_n with initial temperature x_y_z
        rfft = convert.(ComplexF64,rfft) #casting array to ComplexF64
        psd = fourierAnalysis.compute_psd(rfft) #array with the psd THE RFFT[M_n]
        push!(psd_array,psd)

        if i == 1
            str_dashed_temp1 = get_str_dashed_temp(str_rfft_path)
            push!(str_temp_array,str_dashed_temp1)   
        end

        if i == length(str_rfft_path)
            str_dashed_temp2 = get_str_dashed_temp(str_rfft_path)
            push!(str_temp_array,str_dashed_temp2)  
        end    
    end
    
    #string manipulations
    full_file_path *= "psd_$(str_dashed_temp_array[1])_to_$(str_dashed_temp_array[2])_r$(run).pdf" 
    replace_dict = Dict("_" => ".")
    str_temp_array = replace_with_dict.(str_dashed_temp,replace_dict)

    #sampling frecuencies
    f = sampling_freq_arr(str_rfft_path,temp_array)
    
    #= plot styling =#
    plt = plot(f, psd_array, label=L"PSD \ \left( f \right)") #plot reference 
    

    title!("PSD associated to a ts with init temps from $(str_temp_array[1]) to $(str_temp_array[n])")
    xlabel!(L"f")
    ylabel!("power density spectra")

    #= file saving  =#
    savefig(plt, full_file_path)    
end
end #end of module 
