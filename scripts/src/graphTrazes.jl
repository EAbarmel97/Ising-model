module graphTrazes
export save_traze, PlottingException, graph_and_write_over_file!, plot_mean_magn

using Plots
using LaTeXStrings 

include("utilities.jl")
using .utilities: get_array_from_txt, mean_value, neglect_N_first_from_array!

include("exceptions.jl")
using .exceptions: PlottingException

include("ising.jl")
using .ising: CRITICAL_TEMP

#= Function to save the traces of the time series contained in .txt files =#
function save_traze(dir_to_save::AbstractString, 
                        file_path::AbstractString)
    mean = utilities.mean_value(file_path)
    time_series = utilities.get_array_from_txt(file_path)
    x = collect(0:(length(time_series)-1))
    y = time_series
    plt = plot(x, y, label= L"M_n") #plot reference 
    hline!(plt, [mean, mean], label=L"\overline{M}_n",linewidth=3)
    ylims!(-1.0, 1.0)
    xlims!(0, length(time_series))
    xlabel!(L"n")
    ylabel!(L"M_n")
    savefig(plt, dir_to_save) #saving plot reference as a file with pdf extension at a given directory  
end

#= 
Function to write over file and plot the time series contained in each of the all_simulations subdirectories 
meeting a given regex
=#
function graph_and_write_over_file!(dir_names :: AbstractArray, simuls_dir :: AbstractString,
                    file_to_write :: AbstractString, rgx :: Regex)

    curr_dir = pwd()
    GRAPHS_DIR =  joinpath(curr_dir * "graphs")
    GRAPHS_AUTOMATED_DIR = joinpath(curr_dir , "graphs", "automated")

    if !isdir(GRAPHS_AUTOMATED_DIR)
        mkpath(GRAPHS_AUTOMATED_DIR)
    end

    #filtering all file names that match the given regex 
    filtered_array = filter(str -> contains(str, rgx), dir_names)
    
    if isempty(filtered_array)
        throw(exceptions.PlottingException("impossible to graph the given array of temperatures!"))
    end     
    
    for i in eachindex(filtered_array)
        aux_dir_name = filtered_array[i]
        mean_per_temp = 0 
        temp_abs_dir = joinpath(curr_dir,simuls_dir,filtered_array[i],"magnetization") #abs path to the simulations at a given temp 
        num_runs = length(readdir(temp_abs_dir)) #number of runs contained in a given simulations dir
        #creation of the sub dir that will be populated with graphs of magnetization times series at a given temp and run 
        at_temp = replace(filtered_array[i],"simulations_" => "") #stringified temperature 

        if contains(simuls_dir,"/automated/")
            at_temp_dir = joinpath(GRAPHS_AUTOMATED_DIR, at_temp)
            mkpath(at_temp_dir)
        else
            at_temp_dir = joinpath(GRAPHS_DIR, at_temp)
            mkpath(at_temp_dir)   
        end
        
        for run in 1:num_runs
            #.txt magnetization time seris file at a given run
            aux_dir = simuls_dir * filtered_array[i] * "/magnetization/global_magnetization_r$run.txt"
            abs_mean_val = abs(utilities.mean_value(aux_dir))
            mean_per_temp += abs_mean_val/num_runs #mean magnetization at a given temp 
            
            
            #building the graph file name 
            aux_graph_file_name = replace(aux_dir_name,"simulations_T_" => "magnetization_ts_")
            aux_graph_file_name *= "_r$run.pdf"

            if contains(aux_dir,"automated")
                aux_graph_full_name =  joinpath(GRAPHS_AUTOMATED_DIR, at_temp, aux_graph_file_name)
            else
                aux_graph_full_name = joinpath(GRAPHS_DIR, at_temp, aux_graph_file_name)
            end
    
            if isfile(aux_dir) #plot if file exists
                #saves under a common sub dir all simulations with the same temperature
                save_traze(aux_graph_full_name, aux_dir)
            end
        end

        #appending mean value of the global magn time series       
        str_mean_val = "$mean_per_temp"
        aux_temp = replace(aux_dir_name, "simulations_T_" => "", "_" => ".") #getting the temperature
        str_to_append = "$(aux_temp)," * str_mean_val * "\n"
        mean_vals_file = open(file_to_write, "a+")
        write(mean_vals_file, str_to_append) 
        close(mean_vals_file)
    end  
end

#= method to plot custom csv file containing mean magn at its corresponding temp =#
function plot_mean_magn(file_dir :: AbstractString, dir_to_save :: AbstractString)
    temps = []
    mean_magns = []
    # preprocesing custom csv file 
    mean_magn_file = open(file_dir, "r+")
    arr_str = readlines(mean_magn_file) 
    neglect_N_first_from_array!(arr_str,1) #discarting the headers

    for i in eachindex(arr_str)
        substr_temp_and_mean_magn_arr = split(arr_str[i],",")
        stringified_temp = string(substr_temp_and_mean_magn_arr[1])
        stringified_mean_magn = string(substr_temp_and_mean_magn_arr[2])
        temp = utilities.parse_int_float64(Float64,stringified_temp) 
        mean_magn = utilities.parse_int_float64(Float64,stringified_mean_magn)
        push!(temps,temp)
        push!(mean_magns, mean_magn)
    end

    if isfile(file_dir)
        plt = plot(temps, mean_magns, label = L"\overline{M}_n")
        ylims!(0.0, 1.0)
        xlims!(0,3.5)
        vline!(plt, [ising.CRITICAL_TEMP, ising.CRITICAL_TEMP], label=L"T_c", linewidth=1, fillalpha=0.02)
        xlabel!(L"T")
        ylabel!("mean magnetization")
        savefig(plt, dir_to_save) #saving plot reference as a file with pdf extension at a given directory 
    end 
    close(mean_magn_file)
end
end #end of modules