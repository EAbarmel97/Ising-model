module utilities

export binarysearch, quicksort!,swap!, parse_int_float64

include("../src/ising.jl")
using .ising: isingModel,CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS,GLAUBER_DYNAMICS

#=Quick Sort =#

#helper function for swaping values at diferent index locations returning array
function swap!(val1::Int, val2::Int, obj :: Union{Array{Float64,1}, Array{Int,1}})
    temp = obj[val1]
    obj[val1] = obj[val2]
    obj[val2] = temp
    setindex!(obj,obj[val1],val1)
    setindex!(obj,obj[val2],val2)
end

#= partition helper function. Orders the array in a way that all elements to the left of the pivot val are less 
   than the pivot. Elements to the right of the pivot are greather. =#
function partition!(items::Array{Float64,1}, left_pointer::Int, right_pointer::Int)::Int
    pivot_index = right_pointer
    right_pointer -= 1

    while true
        while items[left_pointer] < pivot
            left_pointer += 1
        end

        while items[right_pointer] > pivot
            right_pointer -= 1
        end

        if aux >= right_pointer
            break
        else
            swap!(left_pointer,right_pointer,items)
        end
        swap!(left_pointer,pivot_index,items)
        return left_pointer
    end
end

function qsort(data::Array{Float64,1}, left_pointer::Int, rigth_pointer::Int) :: nothing
    if left_pointer < rigth_pointer
        p = partition(data, left_pointer, rigth_pointer)
        qsort!(data, p, length(data))
        qsort!(data, p + 1, rigth_pointer)
    end
    return nothing 
end

function quicksort!(data::Array{Float64,1})::Array{Doubles64,1}
    qsort(data, 1, length(data))
    return data
end

#=Binary Search =#

#=function to perform binary search. If item not present in array, -1 is returned=#
function bsearch(data::Array{Float64,1}, left_pointer::Int, right_pointer::Int, item ::Float64)::Int
    
    if left_pointer  == right_pointer
        if data[left_pointer] == item
            return left_pointer
        else
            return -1 #item not present en the array in the denerate case 
        end    
    end

    mid_val = (left_pointer + right_pointer)/2
    center = floor(Int,mid_val)
    #item is not in the array 
    if (data[center] < data[left_pointer] )|| (data[right_pointer] < data[center])
        return -1 
    end 
    
    if data[center] == item #item in the center of the array 
        return center
    elseif data[center] > item
        return bsearch(data, left_pointer,center-1, item) #item is in the left subarray 
    else
        return bsearch(data,center+1,right_pointer, item) #item is in the right subarray 
    end    

end

function binarysearch(data::Array{Float64,1}, item :: Float64) :: Int
    return bsearch(data, 1, length(data), item)
end

#= Function to parse an int of float64 at once=#
function parse_int_float64(tp :: Union{Type{Float64},Type{Int}},
    parsing_string :: String) :: Union{Float64,Int}
    try
        parsed_result = Base.parse(tp, parsing_string)
        return parsed_result
    catch e
        isa(e, ArgumentError) 
        printstyled(stderr,"ERROR: imposible to parse a type $tp from '$parsing_string'", 
        bold=true, color=:red) #customized error message 
        println(stderr) 
    end
end

#= TO DO: implement a funtion that burns out the first N simulations =#
#= Given a csv obtained from simulatin, N of first simulations are discarted=#
function prune_N_first_simul(arr :: IOBuffer, N ::Int ) :: nothing
    return nothing
end

end