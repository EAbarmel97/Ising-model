module GraphTrazes
export save_traze, PlottingException, graph_and_write_over_file!, plot_mean_magn

using Plots
using Statistics
using LaTeXStrings 

include("Ising.jl")
using .Ising: CRITICAL_TEMP

include("Exceptions.jl")
using .Exceptions: PlottingException

include("utils/utilities.jl")
using .utilities: get_array_from_txt, mean_value, median_value, neglect_N_first_from_array!,parse_int_float64

include("utils/paths.jl")

#= Function to save the traces of the time series contained in .txt files =#
function save_traze(dir_to_save::String, file_path::String)
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

# Function 1: Create Automated Graphs Directory
function create_automated_graphs_dir()
    curr_dir = pwd()
    GRAPHS_DIR = joinpath(curr_dir, "graphs")
    GRAPHS_AUTOMATED_DIR = joinpath(curr_dir, "graphs", "automated")

    if !isdir(GRAPHS_AUTOMATED_DIR)
        mkpath(GRAPHS_AUTOMATED_DIR)
    end
end

# Function 2: Filter and Check Array
function filter_and_check_array(dir_names::AbstractArray, rgx::Regex)
    filtered_array = filter(str -> contains(str, rgx), dir_names)
    if isempty(filtered_array)
        throw(Exceptions.PlottingException("Impossible to graph the given array of temperatures!"))
    end

    return filtered_array
end

# Function 3: Calculate Median Magnetization
function calculate_median_magnetization(temp_abs_dir::String, num_runs::Int64)
    magnetization_per_run = Float64[]
    for run in 1:num_runs
        aux_dir = joinpath(temp_abs_dir, "global_magnetization_r$run.txt")
        abs_mean_val = abs(utilities.median_value(aux_dir))
        push!(magnetization_per_run, abs_mean_val)
    end

    return median(magnetization_per_run)
end

# Function 4: Create Temperature Directory
function create_temperature_directory(simuls_dir, at_temp, GRAPHS_DIR, GRAPHS_AUTOMATED_DIR)
    if contains(simuls_dir, "/automated/")
        at_temp_dir = joinpath(GRAPHS_AUTOMATED_DIR, at_temp)
        mkpath(at_temp_dir)
    else
        at_temp_dir = joinpath(GRAPHS_DIR, at_temp)
        mkpath(at_temp_dir)
    end
    return at_temp_dir
end

# Function 5: Save Graphs
function save_graphs(temp_abs_dir, aux_dir_name, run, at_temp, GRAPHS_DIR, GRAPHS_AUTOMATED_DIR)
    aux_graph_file_name = replace(aux_dir_name, "simulations_T_" => "magnetization_ts_")
    aux_graph_file_name *= "_r$run.pdf"

    if contains(temp_abs_dir, "automated")
        aux_graph_full_name = joinpath(GRAPHS_AUTOMATED_DIR, at_temp, aux_graph_file_name)
    else
        aux_graph_full_name = joinpath(GRAPHS_DIR, at_temp, aux_graph_file_name)
    end

    if isfile(joinpath(temp_abs_dir, "global_magnetization_r$run.txt"))
        save_traze(aux_graph_full_name, joinpath(temp_abs_dir, "global_magnetization_r$run.txt"))
    end
end

# Function 6: Get Temperatures
function get_temperatures(filtered_array)
    temperatures = Float64[]
    for i in eachindex(filtered_array)
        aux_temp = replace(filtered_array[i], "simulations_T_" => "", "_" => ".")
        temp = utilities.parse_int_float64(Float64, aux_temp)
        push!(temperatures, temp)
    end
    return temperatures
end

# Function 7: Main Function
function graph_and_write_over_file!(dir_names::AbstractArray, simuls_dir::String, file_to_write::String, rgx::Regex)
    create_automated_graphs_dir()
    filtered_array = filter_and_check_array(dir_names, rgx)
    
    temperatures_median_magn = Dict{Float64, String}()
    
    for i in eachindex(filtered_array)
        aux_dir_name = filtered_array[i]
        temp_abs_dir = joinpath(pwd(), simuls_dir, aux_dir_name, "magnetization")
        num_runs = length(readdir(temp_abs_dir))
        at_temp = replace(aux_dir_name, "simulations_" => "")
        at_temp_dir = create_temperature_directory(simuls_dir, at_temp, GRAPHS_DIR, GRAPHS_AUTOMATED_DIR)
        median_per_temp = calculate_median_magnetization(temp_abs_dir, num_runs)
        temperatures_median_magn[temp] = "$median_per_temp"
        
        for run in 1:num_runs
            save_graphs(temp_abs_dir, aux_dir_name, run, at_temp, GRAPHS_DIR, GRAPHS_AUTOMATED_DIR)
        end
    end
    temperatures = get_temperatures(filtered_array)
    return temperatures_median_magn
end


#= method to plot custom csv file containing mean magn at its corresponding temp =#
function plot_mean_magn(file_dir::String, dir_to_save::String)
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
        vline!(plt, [Ising.CRITICAL_TEMP, Ising.CRITICAL_TEMP], label=L"T_c", linewidth=1, fillalpha=0.02)
        xlabel!(L"T")
        ylabel!("mean magnetization")
        savefig(plt, dir_to_save) #saving plot reference as a file with pdf extension at a given directory 
    end 
    close(mean_magn_file)
end
end #end of modules