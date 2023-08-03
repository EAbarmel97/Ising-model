module exceptions
export IlegalChoiceException, PlottingException
#= Exception to be throwned if user makes a bad election  =#
mutable struct IlegalChoiceException <: Exception
    msg::String
end

#showerror function override
Base.showerror(io :: IO, e :: IlegalChoiceException) = print(io, "$(e.msg)")

#= Custom exception to be prompted to the user =#
mutable struct PlottingException <: Exception
    msg::String
end

#showerror function override 
Base.showerror(io :: IO, e :: PlottingException) = print(io, "$(e.msg)")
end #end of module