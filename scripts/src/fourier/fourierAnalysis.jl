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

#= Function to write over a .txt file a vector =#
function write_rfft(arr :: Array{ComplexF64,1}, destination_dir :: AbstractStrings, 
    at_temp :: Float64)
    rounded_temp = round(at_temp, digits=2)
    str_rounded_temp = replace!("$rounded_temp","." => "_")
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
function compute_psd(arr :: Array{ComplexF64,1}) :: Array{Float64,1}
    return abs2.(arr)
end

#= Function to determine the array containing sampling frecuencies =#
function sampling_freq_arr(rfft_paths :: Union{AbstractString, AbstractArray},
    arr_str_temp :: AbstractArray = []) :: Array{Float64,1}
    if typeof(rfft_paths) <: AbstractString  && isempty(arr_str_temp)
        rfft = utilities.get_array_from_txt(rfft_paths)
        return rfftfreq(length(rfft))
    end
    #input is an aray of directories
    if typeof(rfft_paths) <: AbstractArray && !isempty(arr_str_temp)
        aux_rfft_paths = string.(rfft_paths,"/rfft_global_magnetization",arr_str_temp,".txt")
        aux_rfft = utilities.get_array_from_txt(aux_rfft_paths[1]) 
        return rfftfreq(length(aux_rfft))
    end
end

#= 
Function to plot the psd with option to plot it in the same canvas providing a dir where
.txt files with rfft to read from are available or to plot from single text. By default different_canvas is 
set to be true.

NOTE: 
Plot is diplayed in a frecuency vs pwd manner 
=#

function plot_psd(rfft_paths :: Union{AbstractString, AbstractArray}, dir_to_save :: AbstractString,
    different_canvas :: Bool = true)

    graph_file_name = dir_to_save * "psd.pdf"
    #individual psd is plotted
    if  typeof(rfft_paths) <: AbstractString  
        x = sampling_freq_arr(rfft_paths)
        psd_rfft = compute_psd(rfft)
        #= plot styles =#
        plt = plot(x, psd_rfft, label=L"$\lVert DFT[M_n] \rVert^2$") #plot reference 
        
        title!(L"Power density spectrum from $M_n$, initital temperature")
        xlabel!(L"f")
        ylabel!(L"PSD(f)")
        savefig(plt, graph_file_name)
    end

    #multiples psd's are plotted on the same canvas
    if  typeof(rfft_paths) <: AbstractArray 
        psds_array = []
        arr_str_temps = []
        
        for i in eachindex(rfft_paths)
            str_temp = replace(rfft_paths[i], "simulations_" => "")
            push!(arr_str_temps,str_temp)
        end
        
        aux_rfft_paths = string.(rfft_paths,"/rfft_global_magnetization", 
                                arr_str_temps, ".txt") 

        x =  sampling_freq_arr(aux_rfft_paths)
        #multiples psd's are plotted one by one 
        if different_canvas
            for i in eachindex(aux_rfft_paths)
                rfft = utilities.get_array_from_txt(aux_rfft_paths[i])
                psd_rfft = compute_psd(rfft)
                plt = plot(x,psd_rfft, label=L"DFT(X(t))") #plot reference 
                #= plot styles =#
                xlabel!(L"\omega") 
                ylabel!("psd")
                savefig(plt, graph_file_name) 
            end    
        end 

        for i in eachindex(aux_rfft_paths)
            rfft = utilities.get_array_from_txt(rfft_paths[i])
            psd_rfft = compute_psd(rfft)
            push!(psds_array)
        end
    
        plt = plot(x, psds_array) #plot reference 
        #= plot styles =#
        fillalpha!(0.2)
        xlabel!("frec")
        ylabel!("psd")
        savefig(plt, graph_file_name)   
    end    
end
end #end of module 
