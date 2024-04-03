module Exceptions
export IlegalChoiceException, PlottingException, NotIntegerException

struct IlegalChoiceException <: Exception
    msg::String
end

struct PlottingException <: Exception
    msg::String
end
 
struct NotIntegerException <: Exception
    msg::String
end

function Base.showerror(io::IO, e::T) where {T <: Exception}
    print(io, "$(e.msg)")
end

end #end of module