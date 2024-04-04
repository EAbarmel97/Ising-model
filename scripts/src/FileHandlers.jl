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
                throw(ArgumentError("Impossible to parse from $line a type $T"))
            end
        end
        
        return arr  
    catch e
        throw(SystemError("ERROR: There's no such file with path $(text_handler.file_path)"))
    end
end
 
function parse_data(T::Type, text_handler::TextFileHandler)::Dict{Symbol,Vector{T}}
    data_frame = Dict{Symbol, Vector{T}}()
    io = text_handler.file_path
    try 
        open(io, "r") do file
            lines = readlines(file)
            headers = Symbol.(split(lines[1], ","))

            for i in eachindex(lines[2:end])
                fields = split(lines[i], ",")
                for (header, field) in zip(headers, fields)
                    try
                        entry = parse(T, field) 
                        push!(get!(data_frame, header, Vector{T}()), entry)
                    catch e
                        throw(ArgumentError("Unable to parse '$field' as type $T"))
                    end
                end 
            end
        end
        return data_frame
    catch e
        throw(SystemError("ERROR: File not found at path $(text_handler.file_path)"))
    end
end

struct PDFFileHandler <: AbstractFileHandler
    file_path::String
end

function save_file(plt::Plots.Plot{Plots.GRBackend},pdf_handler::PDFFileHandler)
    if !isfile(pdf_handler.file_path)
        savefig(plt, pdf_handler.file_path)
    end

    @error "File $(pdf_handler.file_path) already exists! So, not saving"   
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