module utilities

export binarysearch, quicksort!,swap!, parse_int_float64, prune_N_first

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

#= Given a path to a .txt file containig a time series, N of first simulations are discarted=#
function prune_N_first(full_path_file :: String, N ::Int)    
    time_series = [] #array containg the pruned time series

    file_name = replace(full_path_file,".txt" => "")
    touch("$(file_name)_pruned_$(N)_first_simulations"*".txt")
    
    try
        file = open(full_path_file,"r")
        string_array = readlines(file) # vector with all lines
        for i in 1:eachindex(string_array)
            if i > N 
               push!(time_series,parse_int_float64(Float64,string_array[i])) #parsing stringified i-th observation

               new_file = open(file_name,"a+")
               write(new_file,"$(time_series[i])") #i-th observation of the pruned time series 
               close(new_file)
            end 
        end 
    catch e
        isa(e,SystemError)
        printstyled(stderr,"ERROR: There's no such file with path $(path)", 
        bold=true, color=:red) #customized error message 
        println(stderr)
    end
end
end