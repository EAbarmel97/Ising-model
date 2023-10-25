include("src/IsingModel.jl")
using .IsingModel:do_simulations

include("src/utils/utilities.jl")
using .utilities: get_ARGS

const ARGS_INTERACTIVE = utilities.get_ARGS() #ask user to 

function main()
    IsingModel.do_simulations(ARGS_INTERACTIVE,false) 
end

main()