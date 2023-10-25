#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../src/GraphTrazes.jl")
using .GraphTrazes: save_traze, graph_and_write_over_file!, plot_mean_magn

include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: write_rfft

include("../src/utils/paths.jl")

#auxiliary file directory definitions
const AUTOMATED_SIMULS_DIR_NAMES = readdir(AUTOMATED_SIMULS_DIR)

#= file manipulations =#
touch("mean_magn_automated.txt") #creation of file containing all mean global magnetization
file_to_write_over = joinpath("median_magn_automated.txt")

#= plotting each temperature group and appending it over the file median_magn_automated.txt =#
GraphTrazes.graph_and_write_over_file!(AUTOMATED_SIMULS_DIR_NAMES,AUTOMATED_SIMULS_DIR,file_to_write_over)

#= ploting mean magnetization vs temperature =#
graphTrazes.plot_mean_magn("median_magn_automated.txt", "median_magn_automated.pdf")