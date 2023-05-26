module graphTrazes
export save_traze 

include("utilities.jl")
using .utilities: get_array_from_txt

using Plots

#= Function to save the traces of the time series contained in .txt files =#
function save_traze(dir_to_save :: AbstractString, file_path :: AbstractString)
    time_series = utilities.get_array_from_txt(file_path)
    if !isempty(time_series)
        plt = plot(time_series, xlabel = "observations") #plot reference 
        savefig(plt,dir_to_save) #saving plot reference as a file with pdf extension at a given directory  
    end    
end    
end #end of module