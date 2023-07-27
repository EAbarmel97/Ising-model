module utilities

export swap!, parse_int_float64, get_array_from_txt, mean_value

include("../src/ising.jl")
using .ising: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS

#= Function for swaping values at diferent index locations returning array =#
function swap!(val1::Int, val2::Int, obj::Union{Array{Float64,1},Array{Int,1}})
    temp = obj[val1]
    obj[val1] = obj[val2]
    obj[val2] = temp
    setindex!(obj, obj[val1], val1)
    setindex!(obj, obj[val2], val2)
end

#= Function to parse an int of float64 at once=#
function parse_int_float64(tp::Union{Type{Float64},Type{Int}},
    parsing_string::String)::Union{Float64,Int}
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

#= Function to neglect  the fist N entries from an array =#
function neglect_N_first_from_array!(arr::AbstractArray, first_N::Int)
    try
        arr = arr[(first_N+1):length(arr)]
    catch e
        isa(e, BoundsError)
        printstyled(stderr, "ERROR: cannot neglect first $(first_N)",
            bold=true, color=:red) #customized error message 
        println(stderr)
    end
end

#= Gets an array of strings form a .txt file =#
function get_str_array(file_path::AbstractString)::Array{AbstractString,1}
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
function get_array_from_txt(file_path::AbstractString, prune_first_N=0::Int)::Array{AbstractFloat,1}
    time_series = []

    stringified_array = get_str_array(file_path) #attemps to get an array with the lines of the .txt file
    neglect_N_first_from_array!(stringified_array,prune_first_N) #attemps prunning first N elements from array 
    for i in eachindex(stringified_array)
        push!(time_series, parse_int_float64(Float64, stringified_array[i]))
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
end #end of module
