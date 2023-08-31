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

#= Function to compute the power spectral density =#
function compute_psd(arr :: Array{T,1}) :: Array{Float64,1} where T <: Complex  
    return abs2.(arr)
end

#= Function to determine the array containing sampling frecuencies =#
function sampling_freq_arr(rfft_paths :: Union{AbstractString, AbstractArray},
    arr_str_temp :: AbstractArray = []) :: Array{Float64,1}
    if typeof(rfft_paths) <: AbstractString  && isempty(arr_str_temp)
        rfft = utilities.get_array_from_txt(ComplexF64,rfft_paths)
        return rfftfreq(length(rfft))
    end
    #input is an aray of directories
    if typeof(rfft_paths) <: AbstractArray && !isempty(arr_str_temp)
        #array with full paths to to the stringified rftts files
        aux_rfft_paths = string.(rfft_paths,"/rfft_global_magnetization",arr_str_temp,".txt")
        aux_rfft = utilities.get_array_from_txt(ComplexF64,aux_rfft_paths[1]) 
        return rfftfreq(length(aux_rfft))
    end
end

#= 
Module method for plotting psd wuth options to plot several psd on the same canvas or one by one, providing one generic 
under which all psd will be saved. By default the argument different_canvas is set to be true.

NOTE: full path/s must be given as arguments 
=#
function plot_psd(str_rfft :: Union{AbstractString, AbstractArray}, destination_dir :: AbstractString, different_canvas  = true :: Bool)
    #just one psd will be plotted
    if typeof(str_rfft) <: AbstractString && different_canvas
        #array of sampling frequencies
        f = sampling_freq_arr(str_rfft) 
        rfft = utilities.parse_complex(ComplexF64,str_rfft) #rfft of the M_n with initial temperature x_y_z
        psd = fourierAnalysis.compute_psd(rfft) #array with the psd THE RFFT[M_n]
        
        #getting stringified temperature
        if contains(str_rfft,"all_simualtions/automated/")
            str_temp = replace(str_rfft,r"all_simualtions/automated/simualtions_T_\d_\d_\d" =>"",
                                        "/fourier/rfft_global_magnetization_" => "", r"_\d.txt" => "")
            full_file_path = destination_dir * AUTOMATED_PSD_GRAPHS *"psd_$(str_temp).pdf"
        else
            str_temp = replace(str_rfft,r"all_simualtions/simualtions_T_\d_\d_\d" =>"",
                                        "/fourier/rfft_global_magnetization_" => "", r"_\d.txt" => "")
            full_file_path = destination_dir * PSD_GRAPHS *"psd_$(str_temp).pdf"
        end

        #= plot styling =#
        plt = plot(f, psd, label= L"{\lVert RFFT\left[ M_n \right]\rVert}^2" ) #plot reference 
        title!(L" PSD with initial temperature $%$(str_temp)$")
        xlabel!(L"f")
        ylabel!(L"PSD\left( f \right)")
        #= file saving  =#
        savefig(plt, full_file_path) 
    end    
end
end #end of module 
