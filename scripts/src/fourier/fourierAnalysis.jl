module fourierAnalysis
export  compute_rfft, compute_psd, write_rfft,plot_psd

using FFTW
using Plots
using LaTeXStrings

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

#= Function to determine the array containing sampling frecuencies =#
function sampling_freq_arr(rfft_paths :: Union{AbstractString, AbstractArray},
    arr_str_temp :: AbstractArray = []) :: Array{Float64,1}
    if typeof(rfft_paths) <: AbstractString  && isempty(arr_str_temp)
        rfft = utilities.get_array_from_txt(ComplexF64,rfft_paths)
        return fftfreq(length(rfft))
    end
    #input is an aray of directories
    if typeof(rfft_paths) <: AbstractArray && !isempty(arr_str_temp)
        #array with full paths to to the stringified rftts files
        aux_rfft_paths = string.(rfft_paths,"/rfft_global_magnetization",arr_str_temp,".txt")
        aux_rfft = utilities.get_array_from_txt(ComplexF64,aux_rfft_paths[1]) 
        return fftfreq(length(aux_rfft))
    end
end

#= 
Module method for plotting psd wuth options to plot several psd on the same canvas or one by one, providing one generic 
under which all psd will be saved. By default the argument different_canvas is set to be true.

NOTE: the path/s to the rfft .txt files need to be absulute path/s 
=#
function plot_psd(str_rfft_path :: Union{AbstractString, AbstractArray}, destination_dir :: AbstractString, different_canvas  = true :: Bool)
    curr_dir = pwd()
    full_file_path = destination_dir
    #just one psd will be plotted
    if typeof(str_rfft_path) <: AbstractString && different_canvas
        #array of sampling frequencies
        f = sampling_freq_arr(str_rfft_path)

        rfft = utilities.get_array_from_txt(Complex{Float64},str_rfft_path) #rfft of the M_n with initial temperature x_y_z
        rfft = convert.(ComplexF64,rfft) #casting array to ComplexF64

        psd = fourierAnalysis.compute_psd(rfft) #array with the psd THE RFFT[M_n]
         
        #getting stringified temperature but with "_" as a digit separator
        if contains(str_rfft_path,r"/all_simulations/automated/")
            str_temp = replace(str_rfft_path, curr_dir => "","/all_simulations/automated/simulations_T_" => "",
            r"/fourier/rfft_global_magnetization_\d{1}_\d{1,2}_\d{1}.txt" => "")

        else
            str_temp = replace(str_rfft_path, curr_dir => "","/all_simulations/simulations_T_" => "",
            r"/fourier/rfft_global_magnetization_\d{1}_\d{1,2}_\d{1}.txt" => "")
    
        end
        full_file_path *= "psd_$(str_temp).pdf" 

        #strigified temperature 
        str_temp = replace(str_temp, "_" => ".")

        #= plot styling =#
        plt = plot(f, psd, label=L"PSD \ \left( f \right)") #plot reference 
        title!("PSD associated to a ts with init temp = $(str_temp)")
        xlabel!(L"f")
        ylabel!("power density spectra")
        #= file saving  =#
        savefig(plt, full_file_path)
    end
    
    if typeof(str_rfft_path) <: AbstractArray && different_canvas
        #= TO DO: implement logic when an array of full paths to the stringified rfft are given =#
    end
end
end #end of module 
