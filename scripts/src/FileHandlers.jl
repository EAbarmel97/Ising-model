module FileHandlers
using Plots

include("utils/paths.jl")

export FileHandlerFactory, FileHandler, TextFileHandler
export create_handler, read_data, save_file

abstract type AbstractFileHandler end

struct FileHandler{H <: AbstractFileHandler}
    handler::H
end

struct TextFileHandler <: AbstractFileHandler
    file_path::String
end

function read_data(T::Type, text_handler::TextFileHandler)::Vector{T}
    try 
        lines = String[]
        io = text_handler.file_path
        open(io, "r+") do io
            lines = readlines(io)
        end
       
        arr = T[]
        for line in lines
            try
                parsed_value = parse(T, lines[i])
                push!(arr, parsed_value)
            catch e
                throw(ArgumentError("$line is not of type $T"))
            end
        end
        
        return arr  
    catch e
        throw(SystemError("ERROR: There's no such file with path $(text_handler.file_path)"))
    end
end
 
function parse_data(T::Type, text_handler::TextFileHandler)::Dict{Symbol,Vector{T}}
    data_frame = Dict{:Symbol,Vector{T}}()
    io = text_handler.file_path
    open(io, "r+") do io
        lines = readlines(io)
        headers = Symbol.(lines[1])

        for i in 1:lines[2:end]
            substr_temp_and_mean_magn_arr = split(arr_str[i],",")
            stringified_temp = string(substr_temp_and_mean_magn_arr[1])
            stringified_mean_magn = string(substr_temp_and_mean_magn_arr[2])
            temp = utilities.parse(T,stringified_temp) 
            median_magn = utilities.parse(Float64,stringified_mean_magn)
            push!(temps,temp)
            push!(median_magns, median_magn)
        end
    end    
end

struct PDFFileHandler <: AbstractFileHandler
    file_path::String
end

function save_file(plt::Plots.Plot{Plots.GRBackend},pdf_handler::PDFFileHandler)
    if !isfile(pdf_handler.file_path)
        touch(pdf_handler.file_path)
        savefig(plt, pdf_handler.file_path)
    end

    @error "File $pdf_file already exists! So, not saving"   
end



function create_handler(file_path::String)::FileHandler
    if endswith(file_path, ".txt")
        return FileHandler(TextFile(file_path))
    elseif endswith(file_path, ".pdf")
        return FileHandler(PDFFileHandler(file_path))
    else
        @error "Unsupported type!"
    end
end

end # end of module