include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd, plot_psd

#= auxiliary dir definitions =#
const CURR_DIR = pwd()
const SIMULS_DIR = CURR_DIR * "all_simulations"
const PSD_GRAPHS = "graphs/PSD"

if !isdir(PSD_GRAPHS)
    mkpath(PSD_GRAPHS)
end

DFT = fourierAnalysis.compute_rfft(AUTOMATED_SIMULS_DIR *"simulations_T_1_23/magnetization/global_magnetization_r1.txt")
fourierAnalysis.write_rfft(DFT,AUTOMATED_SIMULS_DIR*"simulations_T_1_23/fourier",1.23)

PSD = fourierAnalysis.compute_psd(DFT)



