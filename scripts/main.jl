include("src/IsingModel.jl")
using .IsingModel:do_simulations

include("src/utilities.jl")
using .utilities: get_ARGS

ARGS = utilities.get_ARGS()

function main()
    IsingModel.do_simulations(ARGS) 
end

main()