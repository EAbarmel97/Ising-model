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
using .utilities: get_array_from_txt, mean_value, median_value, neglect_N_first_from_array!,parse_int_float64,create_graphs_directories,filter_directory_names
using .utilities: count_runs_in_dir

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

# Function 5: Save Graphs
function save_graphs(temp_abs_dir::String, aux_dir_name::String, run::Int64, at_temp)
    aux_graph_file_name = replace(aux_dir_name, "simulations_T_" => "magnetization_ts_")
    aux_graph_file_name *= "_r$run.pdf"

    if contains(temp_abs_dir, "automated")
        aux_graph_full_name = joinpath(AUTOMATED_GRAPHS_DIR, at_temp, aux_graph_file_name)
    else
        aux_graph_full_name = joinpath(GRAPHS_DIR, at_temp, aux_graph_file_name)
    end

    if isfile(joinpath(temp_abs_dir, "global_magnetization_r$run.txt"))
        save_traze(aux_graph_full_name, joinpath(temp_abs_dir, "global_magnetization_r$run.txt"))
    end
end

# Function 6: Get Temperatures
function add_temperature_median_magn_to_dict!(aux_dir_name,temperatures_median_magn)
    
    num_runs = length(readdir(temp_abs_dir))
    
    aux_temp = replace(aux_dir_name, "simulations_T_" => "", "_" => ".")
    temp = utilities.parse_int_float64(Float64, aux_temp)

    median_per_temp = calculate_median_magnetization(temp_abs_dir, num_runs)
    temperatures_median_magn[temp] = "$median_per_temp"
end

function write_header(file_name::String,simuls_dir::String)
    open(file_name,"w+") do io
        if contains(simuls_dir,"/automated/")
            write(io,"temp,median_magn_automated\n") 
        else
            write(io,"temp,median_magn_automated\n") 
        end 
    end   
end 

function write_over_file!(file_to_write::String,temperatures_median_magn::Dict{Float64, String},simuls_dir::String)
    write_header(file_to_write,simuls_dir)

    temperatures = sort!(keys(temperatures_median_magn))

    for i in eachindex(filtered_array)
        temp = temperatures[i]
        if i != length(filtered_array)
            str_to_append = string("$temp,", temperatures_median_magn[temp],"\n")    
        else
            str_to_append = string("$temp,", temperatures_median_magn[temp]) 
        end

        open(file_to_write, "a+") do io 
            write(io, str_to_append) 
        end
    end
end

# Function 7: Main Function
function graph_and_write_over_file!(dir_names::AbstractArray, simuls_dir::AbstractString, file_to_write::AbstractString, rgx::Regex)
    utilities.create_graphs_directories()
    filtered_array = utilities.filter_directory_names(dir_names, rgx)
    
    temperatures_median_magn = Dict{Float64, String}()
    
    for i in eachindex(filtered_array)
        aux_dir_name = filtered_array[i]
        num_runs = count_runs_in_dir(simuls_dir,aux_dir_name)
        at_temp = replace(aux_dir_name, "simulations_" => "")
        utilities.create_temperature_directory(at_temp, simuls_dir)

        add_temperature_median_magn_to_dict!(aux_dir_name,temperatures_median_magn)

        for run in 1:num_runs
            save_graphs(temp_abs_dir, aux_dir_name, run, at_temp)
        end
    end
    
    write_over_file!(file_to_write,temperatures_median_magn,simuls_dir)
end


#= method to plot custom csv file containing mean magn at its corresponding temp =#
function plot_mean_magn(file_dir::String, dir_to_save::String)
    temps = []
    mean_magns = []

    # preprocesing custom csv file 
    open(file_dir, "r+") do io
        arr_str = readlines(io) 
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
    end    
end
end #end of modules