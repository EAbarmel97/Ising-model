include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd, plot_psd

#= auxiliary dir definitions =#
const CURR_DIR = pwd()
const SIMULS_DIR = CURR_DIR * "all_simulations"
const PSD_GRAPHS = "graphs/PSD"

if !isdir(PSD_GRAPHS)
    mkpath(PSD_GRAPHS)
end



