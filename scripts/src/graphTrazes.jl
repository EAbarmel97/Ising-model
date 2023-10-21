module graphTrazes
export save_traze, PlottingException, graph_and_write_over_file!, plot_mean_magn

using Plots
using Statistics
using LaTeXStrings 

include("utilities.jl")
using .utilities: get_array_from_txt, mean_value, median_value, neglect_N_first_from_array!,parse_int_float64

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
function create_graphs_dir_if_not_exists()
    curr_dir = pwd()
    GRAPHS_AUTOMATED_DIR = joinpath(curr_dir, "graphs", "automated")

    if !isdir(GRAPHS_AUTOMATED_DIR)
        mkpath(GRAPHS_AUTOMATED_DIR)
    end
end

function create_temperature_sub_dir(sub_dir:: String)
    #creation of the sub dir that will be populated with graphs of magnetization times series at a given temp and run 
    at_temp = replace(sub_dir,"simulations_" => "")
    
    if contains(simuls_dir,"/automated/")
        at_temp_dir = joinpath(GRAPHS_AUTOMATED_DIR, at_temp)
        mkpath(at_temp_dir)
    else
        at_temp_dir = joinpath(GRAPHS_DIR, at_temp)
        mkpath(at_temp_dir)   
    end
end

function count_runs_in_dir(simuls_dir::Array{String,1})::Int64
    temp_abs_dir = joinpath(curr_dir,simuls_dir,aux_dir_name,"magnetization") #abs path to the simulations at a given temp 
    return length(readdir(temp_abs_dir)) #number of runs contained in a given simulations dir
end

function filter_by_rgx(dirs_to_filter::Array{String,1},rgx::Regex)::Array{String,1}
    filtered_array = filter(str -> contains(str, rgx), dirs_to_filter)
    if isempty(filtered_array)
        throw(exceptions.PlottingException("impossible to graph the given array of temperatures!"))
    end
    return filtered_array 
end

function create_ordered_temperatures_median_magn_dict(dir_names :: AbstractArray)::Dict{Float64,String}
    temperatures = Float64[]
    str_median_magn = String[]
    ordered_dict = Dict{Float64,String}()

    for i in eachindex(dir_names)
        aux_temp = replace(dir_names[i], "simulations_T_" => "", "_" => ".") #getting the temperature
        temp = utilities.parse_int_float64(Float64,aux_temp)
        push!(temperatures,temp)

        median_per_temp = median(magnetization_per_run)
        push!(str_median_magn,"$median_per_temp")
    end

    sort!(temperatures)
    for i in eachindex(dir_names)
        ordered_dict[temperatures[i]] = str_median_magn[i]
    end

    return ordered_dict
end

function write_over_file!(dir_names :: AbstractArray, file_to_write :: AbstractString)
    for i in eachindex(dir_names)
        magnetization_per_run = Float64[]

        num_runs = count_runs_in_dir(simuls_dir)
        for run in 1:num_runs
            #.txt magnetization time seris file at a given run
            aux_dir = joinpath(temp_abs_dir,"global_magnetization_r$run.txt")
            abs_median_val = abs(utilities.median_value(aux_dir))
            push!(magnetization_per_run,abs_median_val)  
        end   
    end
    
    temperature_median_magn_dict = create_ordered_temperatures_median_magn_dict(filtered_arr)

    #appending median value of the global magn time series       
    for key in keys(temp_median_magn_dict)
        open(file_to_write, "a+") do io 
            str_to_append = "$key," * temperature_median_magn_dict[key] * "\n"
            write(io, str_to_append) 
        end
    end            
end

function graph_traze(dir_names, simuls_dir)
    num_runs = count_runs_in_dir(dir_names)
    for i in eachindex(dir_names)
        for i in 1:num_runs
            
        end
        #building the graph file name
        aux_dir = joinpath(temp_abs_dir,"global_magnetization_r$run.txt")
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
end


function graph_and_write_over_file!(dir_names :: AbstractArray, simuls_dir :: AbstractString, file_to_write :: AbstractString)
    rgx_arr = [r"T_0_\d{1,2}", r"T_1_\d{1,2}", r"T_2_\d{1,}", r"T_3_\d{1,2}"]

    create_graphs_dir_if_not_exists()

    try
        for i in eachindex(rgx_arr)
            filtered_arr = filter_by_rgx(dir_names,rgx) 
            write_overfile!(filtered_arr,file_to_write)
            graph_traze(filtered_arr,simuls_dir)
        end   
    catch

    end
end =#

#= 
Function to write over file and plot the time series contained in each of the all_simulations subdirectories 
meeting a given regex
=#

function graph_and_write_over_file!(dir_names :: AbstractArray, simuls_dir :: AbstractString,
                    file_to_write :: AbstractString, rgx :: Regex)

    curr_dir = pwd()
    GRAPHS_DIR =  joinpath(curr_dir, "graphs")
    GRAPHS_AUTOMATED_DIR = joinpath(curr_dir, "graphs", "automated")

    if !isdir(GRAPHS_AUTOMATED_DIR)
        mkpath(GRAPHS_AUTOMATED_DIR)
    end

    #filtering all file names that match the given regex 
    filtered_array = filter(str -> contains(str, rgx), dir_names)
    if isempty(filtered_array)
        throw(exceptions.PlottingException("impossible to graph the given array of temperatures!"))
    end 
    
    temperatures = Float64[]
    temperatures_median_magn = Dict{Float64,String}()

    for i in eachindex(filtered_array)
        median_per_temp = 0 
        magnetization_per_run = Float64[]
        
        aux_dir_name = filtered_array[i]
        temp_abs_dir = joinpath(curr_dir,simuls_dir,aux_dir_name,"magnetization") #abs path to the simulations at a given temp 
        num_runs = length(readdir(temp_abs_dir)) #number of runs contained in a given simulations dir
        #creation of the sub dir that will be populated with graphs of magnetization times series at a given temp and run 

        at_temp = replace(filtered_array[i],"simulations_" => "") #stringified dashed temperature with T_ prefixed

        if contains(simuls_dir,"/automated/")
            at_temp_dir = joinpath(GRAPHS_AUTOMATED_DIR, at_temp)
            mkpath(at_temp_dir)
        else
            at_temp_dir = joinpath(GRAPHS_DIR, at_temp)
            mkpath(at_temp_dir)   
        end
        
        for run in 1:num_runs
            #.txt magnetization time seris file at a given run
            aux_dir = joinpath(temp_abs_dir,"global_magnetization_r$run.txt")
            abs_mean_val = abs(utilities.median_value(aux_dir))
            push!(magnetization_per_run,abs_mean_val)  
        end
    
        for run in 1:num_runs
            #building the graph file name
            aux_dir = joinpath(temp_abs_dir,"global_magnetization_r$run.txt")
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

        aux_temp = replace(aux_dir_name, "simulations_T_" => "", "_" => ".") #getting the temperature
        temp = utilities.parse_int_float64(Float64,aux_temp)
        push!(temperatures,temp)

        median_per_temp = median(magnetization_per_run)

        temperatures_median_magn[temp] = "$median_per_temp"    
    end

    sort!(temperatures)

    #appending median value of the global magn time series       
    for i in eachindex(filtered_array)
        temp = temperatures[i]
        open(file_to_write, "a+") do io 
            str_to_append = "$temp," * temperatures_median_magn[temp] * "\n"
            write(io, str_to_append) 
        end
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