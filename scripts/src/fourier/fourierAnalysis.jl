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
function write_rfft(arr :: Array{ComplexF64,1}, destination_dir :: AbstractStrings)
    file_name = "$destination_dir/rfft_global_magnetization.txt"
    touch(file_name)
    fourier_transform_file = open(file_name,"a+")
    for i in eachindex(arr)
        write(fourier_transform_file,"$(arr[i])")
    end
    close(fourier_transform_file)    
end

#= Function to compute the power spectral density =#
function compute_psd(arr :: Array{ComplexF64,1}) :: Array{Float64,1}
    return abs2.(arr)
end

#= 
Function to plot the psd with option to plot it in the same canvas providing a dir where
.txt files with rfft to read from are available or to plot from single text. By default same_canvas is 
set to be false. 

NOTE: 
Plot is diplayed in a frecuency vs pwd manner 
=#
function plot_psd(rfft_path :: Union{AbstractString, AbstractArray}, dir_to_save :: AbstractString,
    same_canvas :: Bool = false)
    #array containing frecuencies
    x = rfftfreq(length(rfft))
    #individual psd is plotted
    if !same_canvas  && typeof(rfft_path) <: AbstractString  
        rfft = utilities.get_array_from_txt(rfft_path)
        psd_rfft = compute_psd(rfft)
        plt = plot(x, psd_rfft) #plot reference 
        xlabel!("frec")
        ylabel!("psd")
        savefig(plt, dir_to_save)
    end

    #multiples psd's are plotted on the same canvas
    if same_canvas && typeof(rfft_path) <: AbstractArray 
        psds_array = []
        #array containing the full paths to each rfft 
        rfft_path = string.(rfft_path,"/rfft_global_magnetization.txt") 
        for i in eachindex(rfft_path)
            rfft = utilities.get_array_from_txt(rfft_path[i])
            psd_rfft = compute_psd(rfft)
            push!(psds_array)
        end

        plt = plot(x,psds_array) #plot reference 
        #= plot styles =#
        fillalpha!(0.2)
        xlabel!("frec")
        ylabel!("psd ")
        savefig(plt, dir_to_save)                 
    end    
end
end #end of module 