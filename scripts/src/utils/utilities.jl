module utilities
export swap!, parse_int_float64, get_array_from_txt, mean_value, push_arith_progression!,use_temperature_array, TEMPERATURE_INTERVALS, replace_with_dict
export median_value, get_ARGS
export create_simulation_sub_dir,create_fourier_dir,create_simulations_dir_if_not_exist,create_dir_if_not_exists,create_graphs_directories,filter_directory_names

export count_runs_in_dir, count_number_of_directories_maching_rgx
export determines_noise_or_movement,graph_file_path

using  Statistics

include("../Ising.jl")
using .Ising: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS

include("../Exceptions.jl")
using .Exceptions: NotIntegerException

include("paths.jl")


function swap!(val1::Int, val2::Int, obj::Union{Array{Float64,1},Array{Int,1}})
    temp = obj[val1]
    obj[val1] = obj[val2]
    obj[val2] = temp
    setindex!(obj, obj[val1], val1)
    setindex!(obj, obj[val2], val2)
end

#= Function to parse an int of float64 at once=#
function parse_int_float64(tp::Union{Type{Float64},Type{Int}},parsing_string::String)::Union{Float64,Int}
    try
        parsed_result = Base.parse(tp, parsing_string)
        return parsed_result
    catch e
        isa(e, ArgumentError)
        printstyled(stderr, "ERROR: imposible to parse a type $tp from '$parsing_string'",
            bold=true, color=:red) #customized error message 
        println(stderr)
    end
end

#= Function to parse complex input ina string =#
function parse_complex(tp::Type{T}, parsing_string::String) where T <: Complex
    try
        parsed_result = Base.parse(tp, parsing_string)
        return parsed_result
    catch e
        isa(e, ArgumentError)
        printstyled(stderr, "ERROR: imposible to parse a type $tp from '$parsing_string'",
            bold=true, color=:red) #customized error message 
        println(stderr)
    end
end

"""
    neglect_N_first_from_array!(arr::AbstractArray, first_N::Int)::Nothing

Dsicards the fist N entries from an array
"""
function neglect_N_first_from_array!(arr::AbstractArray, first_N::Int)::Nothing
    try
        if first_N >= 1            
            deleteat!(arr, 1:first_N)
        end
    catch e
        isa(e, BoundsError)
        printstyled(stderr, "ERROR: cannot neglect first $(first_N)",
            bold=true, color=:red) #customized error message 
            println(stderr)
    end

    return nothing
end

    #= Gets an array of strings form a .txt file =#
function get_str_array(file_path::AbstractString)::Array{String,1}
    stringified_array = []
    try
        opened_file = open(file_path, "r+")
        stringified_array = readlines(opened_file)
        return stringified_array
    catch e
        isa(e, SystemError)
        printstyled(stderr, "ERROR: There's no such file with path $(file_path)",
            bold=true, color=:red) #customized error message 
        println(stderr)
    end
end

"""
    get_array_from_txt(file_path::AbstractString, prune_first_N=0::Int):: Array{Float64,1}

Gets an array of floats from a .txt file, given its path
"""
function get_array_from_txt(file_path::AbstractString, prune_first_N=0::Int):: Array{Float64,1}
    arr = Float64[]
    
    stringified_array = get_str_array(file_path) #attemps to get an array with the lines of the .txt file
    neglect_N_first_from_array!(stringified_array,prune_first_N) #attemps prunning first N elements from array 
    for i in eachindex(stringified_array)
        push!(arr, parse_int_float64(Float64, stringified_array[i]))
    end
    
    return arr 
end

#= Method to parse an array with elements Float64 or ComplexF64 from a .txt file =#
function get_array_from_txt(tp::Union{Type{Float64}, Type{Complex{Float64}}},file_path::AbstractString, prune_first_N=0::Int)::Array{Union{ComplexF64,Float64},1}

    arr = Union{ComplexF64,Float64}[]
    stringified_array = get_str_array(file_path) #attemps to get an array with the lines of the .txt file
    neglect_N_first_from_array!(stringified_array,prune_first_N) #attemps prunning first N elements from array 

    for i in eachindex(stringified_array)
        if tp == Float64
            push!(arr, parse_int_float64(Float64, stringified_array[i]))
        end 
        
        if  tp == Complex{Float64}
            push!(arr, parse_complex(Complex{Float64}, stringified_array[i]))
        end
    end

    return arr
    
end

# Function 6: Get Temperatures
function get_sorted_temperatures(filtered_array)
    temperatures = Float64[]
    
    for i in eachindex(filtered_array)
        aux_temp = replace(filtered_array[i], "simulations_T_" => "", "_" => ".")
        temp = utilities.parse_int_float64(Float64, aux_temp)
        push!(temperatures, temp)
    end
    
    return sort!(temperatures)
end

#= Function to get the arithmetic mean value of a given time series=#
function mean_value(file_path::AbstractString, prune_first_N=0::Int):: Float64
    sum = 0
    time_series = get_array_from_txt(file_path,prune_first_N)
    for i in eachindex(time_series)
        sum += time_series[i]
    end
    return sum /= length(time_series)
end

function median_value(file_path::AbstractString, prune_first_N=0::Int):: Float64
    time_series = get_array_from_txt(file_path,prune_first_N)
    return median(time_series)
end

#= Function to push fill in values in arithmetic progression on an array given the end points=#
function push_arith_progression!(Ti::Float64, Tf::Float64, delta::Float64, arr::AbstractArray)
    val = cld(Tf-Ti,delta)
    umbral = round(abs(Tf - (Ti + val*delta)), digits=2)
    #= If the umbral lies in the interval (0,0.1], then Tf is taken to be in arithmetic 
       progression with Ti =#
    if umbral <= 0.1 
        for i in 1:val
            new_val = Ti + i*delta
            push!(arr,new_val) 
        end  
    else
        throw(exceptions.NotIntegerException("Illegal choice. $Tf is not in arithmetic progression with respect to $Ti"))
    end  
end

function get_default_temperature_array()::Array{Float64}
    default_array = [0.0]
    push_arith_progression!(0.0,1.0,0.1,default_array)
    push_arith_progression!(1.0,2.2,0.1,default_array)
    push_arith_progression!(2.21,2.26,0.01,default_array)
    push!(default_array,Ising.CRITICAL_TEMP)
    push_arith_progression!(2.27,2.5,0.01,default_array)
    push_arith_progression!(2.6,3.5,0.01,default_array)

    return default_array
end

const TEMPERATURE_INTERVALS = get_default_temperature_array()

#= method to replace a function with a replacement dictionary =#
function replace_with_dict(str :: String, replace_dict :: Dict{T, String} where T <: Any) :: String 
    return replace(str,replace_dict...)
end

function separate_by_dashes(str::String)::String
    return replace(str, "." => "_")
end

function append_2_element_array_or_throw(arr1::AbstractArray,arr2)
    if arr1[1] < arr1[2]
        throw(exceptions.IlegalChoiceException("Ilegal choice. Tf < Ti"))  
    end
    append!(arr2,arr1)    
end

#= Auxiliary functions for user input processing =#

function use_temperature_array() :: Bool
    print("use temperature array?: [y/n]\t")
    usr_input = lowercase(readline())
    
    while true
        if usr_input == "y"
            return true
        end
        
        if usr_input == "n"
            return false
        end
        
        usr_input = lowercase(readline())
    end
end

function get_ARGS()::Array{Union{Int64,Float64},1}
    ARGS = Union{Int64,Float64}[]
    
    println("Provide initial and final temperatures. Ti < Tf")
    
    println()
     
    println("Initial temperature:")
    init_temp = parse_int_float64(Float64, readline())
    
    println()
    
    println("Final temperature:")
    final_temp = parse_int_float64(Float64, readline())
    
    append_2_element_array_or_throw([init_temp,final_temp],ARGS)

    println()

    println("Number of runs:")
    num_runs = parse_int_float64(Int, readline())
    push!(ARGS,num_runs)

    println()
     
    println("Increment: ")
    increment = parse_int_float64(Float64, readline())
    push!(ARGS,increment)

    println("Grid size")
    n_grid = parse_int_float64(Int, readline())
    push!(ARGS,n_grid)

    println()
     
    println("Number of generations:")
    num_generations = parse_int_float64(Int, readline())
    push!(ARGS,num_generations)

    return ARGS
end

#= IO handling auxiliary functions =#
function create_automated_simulations_dir_if_not_exists()
    if !isdir(AUTOMATED_SIMULS_DIR)
        mkpath(AUTOMATED_SIMULS_DIR)
    end
end

function create_simulations_dir_if_not_exists()
    if !isdir(SIMULS_DIR)
        mkpath(SIMULS_DIR)
    end
end

function create_simulations_dirs_if_not_exit(is_automated::Bool)
    if is_automated
        create_automated_simulations_dir_if_not_exists()
    else
        create_simulations_dir_if_not_exists()
    end
end

function create_simulation_sub_dir(temp::Float64,is_automated::Bool)::String
    ROUNDED_TEMP = round(temp, digits=2)
    str_temp = replace("$(ROUNDED_TEMP)", "." => "_") #stringified temperature with "." replaced by "_"
    simulations_dir =  abspath(string("simulations_T_", str_temp)) #folder containing simulations all temp str_temp 
    mkpath(simulations_dir)

    return simulations_dir
end

function create_graphs_dir_if_not_exits()
    if !isdir(GRAPHS_DIR_SIMULS)
        mkpath(GRAPHS_DIR_SIMULS)
    end
end

function create_automated_graphs_dir_if_not_exists()
    if !isdir(AUTOMATED_GRAPHS_DIR_SIMULS)
        mkpath(AUTOMATED_GRAPHS_DIR_SIMULS)
    end
end

# Function 1: Create Graphs Directories
function create_graphs_directories(simuls_dir::String)
    if contains(simuls_dir,"automated")
        create_automated_graphs_dir_if_not_exists()
    else       
        create_graphs_dir_if_not_exits()
    end       
end

# Function 4: Create Temperature Directory
function create_temperature_directory(at_temp::String, simuls_dir::String)
    if contains(simuls_dir, "automated")
        #= at_temp_dir = joinpath(AUTOMATED_GRAPHS_DIR, at_temp) =#
        mkpath(joinpath(AUTOMATED_GRAPHS_DIR_SIMULS, at_temp))
    else
        #= at_temp_dir = joinpath(GRAPHS_DIR, at_temp) =#
        mkpath(joinpath(GRAPHS_DIR_SIMULS, at_temp))
    end
end

# Function 2: Filter Directory Names
function filter_directory_names(dir_names, rgx)
    filtered_array = filter(str -> contains(str, rgx), dir_names)
    if isempty(filtered_array)
        throw(Exceptions.PlottingException("Impossible to graph the given array of temperatures!"))
    end
    
    return filtered_array
end

function create_magnetization_dir(simulation_dir::String)::String
    global_magnetization_aux_dir = joinpath(simulation_dir, "magnetization")
    mkpath(global_magnetization_aux_dir)

    return global_magnetization_aux_dir
end

function create_correlated_noise_dir()
    if !isdir(CORRELATED_NOISE_DIR)
        mkpath(CORRELATED_NOISE_DIR)
    end
end

function create_magnetization_time_series_file(file_name, dir_to_save)
    #= Creation of generic .txt files containing global magnetization time series =#
    magnetization_file_path = joinpath(dir_to_save, file_name)
    touch(magnetization_file_name)

    return magnetization_file_path
end

function create_fourier_dir_if_not_exists(simulations_dir::String)::String
    FOURIER_AUTOMATED_DIR = joinpath(simulations_dir, "fourier")
    mkpath(FOURIER_AUTOMATED_DIR)  
end

function create_dir_not_exists(dit_to_create::String)::Nothing
    if !isdir(dir_to_create)
        mkpath(dir_to_create)
    end
    
    return nothing
end    

function count_runs_in_dir(simuls_dir::String, aux_dir_name::String)::Int64
    temp_abs_dir = joinpath(simuls_dir,aux_dir_name,"magnetization") #abs path to the simulations at a given temp 
    return length(readdir(temp_abs_dir)) #number of runs contained in a given simulations dir
end

function count_number_of_directories_maching_rgx(dir_names::Array{String},rgx::Regex)
    return length(filter_directory_names(dir_names,rgx))
end

"""
    psd_graph_file_path(destination_dir::String,beta0::Float64,beta1::Float64)::String

Outputs the file path where the psd associated of the simulated time series will be saved 
"""
function graph_file_path(destination_dir::String,A::Float64,beta::Float64)::String
    full_file_path = joinpath(destination_dir,"ts_log_A_$(separate_by_dashes(string(round(log10(A),digits=2))))_beta_$(separate_by_dashes(string(round(beta,digits=2)))).pdf")

    return full_file_path
end

"""
    psd_graph_file_path(destination_dir::String,beta0::Float64,beta1::Float64)::String

Outputs the file path where the psd associated of the simulated time series will be saved 
"""
function psd_graph_file_path(destination_dir::String,A::Float64,beta::Float64)::String

    full_file_path = joinpath(destination_dir,"psd_log_A_$(separate_by_dashes(string(round(log10(A),digits=2))))_beta_$(separate_by_dashes(string(round(beta,digits=2)))).pdf")

    return full_file_path
end
end #end of module
