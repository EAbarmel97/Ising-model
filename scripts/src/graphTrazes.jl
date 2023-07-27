module graphTrazes
export save_traze

include("utilities.jl")
using .utilities: get_array_from_txt, PlottingException, graph_and_write_over_file

using Plots

#= Custom exception to be prompted to the user =#
struct PlottingException <: Exception
    msg::String
end

#showerror fucntion override
Base.showerror(io :: IO, e :: PlottingException) = print($(e.msg))

#= Function to save the traces of the time series contained in .txt files =#
function save_traze(dir_to_save::AbstractString, file_path::AbstractString)
    mean = utilities.mean_value(file_path)
    time_series = utilities.get_array_from_txt(file_path)
    x = collect(0:(length(time_series)-1))
    y = time_series
    plt = plot(x, y, label="Mg_n") #plot reference 
    hline!([mean, mean])
    ylims!(0.0, 1.0)
    xlims!(0, length(time_series))
    xlabel!("obs")
    ylabel!("magn time series")
    savefig(plt, dir_to_save) #saving plot reference as a file with pdf extension at a given directory  
end

#= Function to write over file and plot the time series contained in each of the all_simulations subdirectories=#
function graph_and_write_over_file!(dir_names :: Array{AbstractString,1}, file_to_write :: AbstractString, rgx :: Regex)
    #filtering all file names that match the given regex 
    filtered_array = filter( str -> contains( str -> rgx), dir_names)
    
    if isempty(filtered_array)
        throw(PlottingException("impossible to graph the given array of temperatures!"))
    end     
    
    for i in eachindex(filtered_array)
        aux_dir =  automated_simuls * "$(filtered_array[i])" * "/magnetization/global_magnetization_r1.txt"
        mean_val = utilities.mean_value(aux_dir) #mean value of the global magnetization time series
        str_mean_val = "$(mean_val)"
        
        #appending mean value of the global magn time series
        aux_temp = replace(filtered_array[i], "simulations_T_" => "T.") #getting the temperature
        str_to_append = "$(aux_temp)," * str_mean_val 
        mean_vals_file = open(file_to_write, "a+")
        write(mean_vals_file, str_to_append) 
        close(mean_vals_file)
        
        aux_graph_file_name = replace(filtered_array[i],".txt" => ".pdf")
        aux_graph_name = curr_dir * "/graphs/" * aux_graph_file_name

        if isfile(aux_dir) #plot if file exists
            save_traze(aux_graph_name, aux_dir)
        end
    end   
end
end #end of modules