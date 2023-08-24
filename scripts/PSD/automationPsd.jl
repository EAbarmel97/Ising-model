include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd, plot_psd

#= auxiliary dir definitions =#
const CURR_DIR = pwd()
const AUTOMATED_SIMULS_DIR = CURR_DIR * "/all_simulations/automated/"
const AUTOMATED_PSD_GRAPHS = "graphs/automated/PSD/"

if !isdir(AUTOMATED_PSD_GRAPHS)
    mkpath(AUTOMATED_PSD_GRAPHS)
end

DFT = fourierAnalysis.compute_rfft(AUTOMATED_SIMULS_DIR *"simulations_T_1_23/magnetization/global_magnetization_r1.txt")
fourierAnalysis.write_rfft(DFT,AUTOMATED_SIMULS_DIR*"simulations_T_1_23/fourier",1.23)

PSD = fourierAnalysis.compute_psd(DFT)
