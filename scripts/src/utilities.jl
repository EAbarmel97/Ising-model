module utilities
export swap!, parse_int_float64, get_array_from_txt, mean_value, push_arith_progression!,use_temperature_array, TEMPERATURE_INTERVALS, replace_with_dict
include("../src/ising.jl")
using .ising: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS

include("../src/exceptions.jl")
using .exceptions: NotIntegerException


#= Function for swaping values at diferent index locations returning array =#
function swap!(val1::Int, val2::Int, obj::Union{Array{Float64,1},Array{Int,1}})
    temp = obj[val1]
    obj[val1] = obj[val2]
    obj[val2] = temp
    setindex!(obj, obj[val1], val1)
    setindex!(obj, obj[val2], val2)
end

#= Function to parse an int of float64 at once=#
function parse_int_float64(tp :: Union{Type{Float64},Type{Int}},
    parsing_string :: String) :: Union{Float64,Int}
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
function parse_complex(tp :: Type{T}, parsing_string :: String) where T <: Complex
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

#= Function to match user input =#
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

#= Function to neglect  the fist N entries from an array =#
function neglect_N_first_from_array!(arr::AbstractArray, first_N::Int)
    try
        deleteat!(arr, 1:first_N)
    catch e
        isa(e, BoundsError)
        printstyled(stderr, "ERROR: cannot neglect first $(first_N)",
            bold=true, color=:red) #customized error message 
            println(stderr)
    end
end

    #= Gets an array of strings form a .txt file =#
function get_str_array(file_path::AbstractString):: Array{String,1}
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

#= Function that gets an array of floats from a .txt file.

NOTE: 
If the .txt file is not on the same hierarchy its absulute path must be provided =#
function get_array_from_txt(file_path::AbstractString, prune_first_N=0::Int):: Array{Float64,1}
    time_series = []
    
    stringified_array = get_str_array(file_path) #attemps to get an array with the lines of the .txt file
    stringified_array = neglect_N_first_from_array!(stringified_array,prune_first_N) #attemps prunning first N elements from array 
    for i in eachindex(stringified_array)
        push!(time_series, parse_int_float64(Float64, stringified_array[i]))
    end
    
    return time_series
end

#= Method to parse an array with elements Float64 or ComplexF64 from a .txt file =#
function get_array_from_txt(tp :: Union{Type{Float64}, Type{Complex{Float64}}},file_path :: AbstractString, 
                prune_first_N=0 :: Int) :: Array{Union{ComplexF64,Float64},1}

    time_series = []
    stringified_array = get_str_array(file_path) #attemps to get an array with the lines of the .txt file
    stringified_array = neglect_N_first_from_array!(stringified_array,prune_first_N) #attemps prunning first N elements from array 

    for i in eachindex(stringified_array)
        if tp == Float64
            push!(time_series, parse_int_float64(Float64, stringified_array[i]))
        end 
        
        if  tp == Complex{Float64}
            push!(time_series, parse_complex(Complex{Float64}, stringified_array[i]))
        end
    end

    return time_series
    
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

#= Function to push fill in values in arithmetic progression on an array given the end points=#
function push_arith_progression!(Ti :: Float64, Tf :: Float64, delta :: Float64,
    arr :: AbstractArray)
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

aux_temps_intervals = [0.0]
push_arith_progression!(0.0,1.0,0.1,aux_temps_intervals)
push_arith_progression!(1.0,2.2,0.1,aux_temps_intervals)
push_arith_progression!(2.2,2.5,0.01,aux_temps_intervals)
push_arith_progression!(2.5,3.5,0.1,aux_temps_intervals)

const TEMPERATURE_INTERVALS = aux_temps_intervals

#= method to replace a function with a replacement dictionary =#
function replace_with_dict(str :: String, replace_dict :: Dict{T, String} where T <: Any) :: String 
    return replace(str,replace_dict...)
end
end #end of module
