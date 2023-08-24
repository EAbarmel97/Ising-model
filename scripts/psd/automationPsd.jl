include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd, plot_psd

#= auxiliary dir definitions =#
const CURR_DIR = pwd()
const AUTOMATED_SIMULS_DIR = CURR_DIR * "/all_simulations/automated/"
const AUTOMATED_PSD_GRAPHS = "graphs/automated/psd/"

if !isdir(AUTOMATED_PSD_GRAPHS)
    mkpath(AUTOMATED_PSD_GRAPHS)
end

DFT = fourierAnalysis.compute_rfft(AUTOMATED_SIMULS_DIR *"simulations_T_2_31/magnetization/global_magnetization_r1.txt")
fourierAnalysis.write_rfft(DFT,AUTOMATED_SIMULS_DIR*"simulations_T_2_31/fourier",2.31)

PSD = fourierAnalysis.compute_psd(DFT)
