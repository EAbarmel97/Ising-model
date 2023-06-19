module graphTrazes
export save_traze

include("utilities.jl")
using .utilities: get_array_from_txt

using Plots

#= Function to get the mean value of a given time series=#
function mean_val(file_path :: AbstractString) :: Float64
    sum = 0
    time_series = utilities.get_array_from_txt(file_path)
    for i in eachindex(time_series)
        sum += time_series[i]
    end
    return sum /= length(time_series)     
end    

#= Function to save the traces of the time series contained in .txt files =#
function save_traze(dir_to_save :: AbstractString, file_path ::AbstractString)
    time_series = utilities.get_array_from_txt(file_path)
    plt = plot(y = time_series, xlabel = "observations", 
    ylabel = "magn time series", ylims = [0.0, 1.0]) #plot reference 
    savefig(plt, dir_to_save) #saving plot reference as a file with pdf extension at a given directory  
end
end #end of module