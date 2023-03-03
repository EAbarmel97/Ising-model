#partition helper function
function partition(items::Array{Double64,1},leftPointer::Int, RightPointer :: Int)
    
end

#quicksort
function quickSort!(data::Array{Double64,1}, leftPointer :: Int, rigthPointer :: Int )
    if leftPointer < rigthPointer
        p = partition(data,leftPointer,pivot)
        quickSort!(data, p, length(data))
        quickSort!(data,p+1,rigthPointer)
    end    
end 

#binary serach
function binarySearch!(params:: Array{})
    return true 
end