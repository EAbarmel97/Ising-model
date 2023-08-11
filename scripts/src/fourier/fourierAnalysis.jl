module fourierAnalysis
export  compute_rfft, compute_psd, write_rfft

using FFTW

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

#= Function to compute the spectral density =#
function compute_psd(arr :: Array{ComplexF64,1})
    return abs2.(arr)
end
end #end of module 