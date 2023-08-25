module fourierAnalysis
export  compute_rfft, compute_psd, write_rfft

using FFTW
using Plots
using LaTeXStrings

include("../utilities.jl")
using .utilities: get_array_from_txt

#= Function wrapper to the FFTW's rfft method =#
function compute_rfft(file_path :: AbstractString) :: Array{ComplexF64,1}
    time_series = utilities.get_array_from_txt(file_path)
    return rfft(time_series) 
end

#= Function to write over a .txt file a vector with the rfft of a signal(time series) =#
function write_rfft(arr :: Array{ComplexF64,1}, destination_dir :: AbstractString, 
    at_temp :: Float64)
    rounded_temp = round(at_temp, digits=2)
    str_rounded_temp = replace("$rounded_temp","." => "_")
    file_name = "$destination_dir/rfft_global_magnetization_$(str_rounded_temp).txt"
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
end #end of module 
