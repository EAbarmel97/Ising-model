include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd 

include("../src/utilities.jl")
using .utilities: get_array_from_txt, parse_int_float64

#= auxiliary dir definitions =#
const CURR_DIR = pwd()
const SIMULS_DIR = CURR_DIR * "all_simulations"
const PSD_GRAPHS = "/graphs/psd/"
const ALL_FOURIER_DIRS = string.(SIMULS_DIR,readdir(SIMULS_DIR),"/fourier/")

if !isdir(PSD_GRAPHS)
    mkpath(PSD_GRAPHS)
end
   


