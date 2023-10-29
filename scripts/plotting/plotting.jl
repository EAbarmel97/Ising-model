#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../src/GraphTrazes.jl")
using .GraphTrazes: save_traze, graph_and_write_over_file!, plot_median_magn

include("../src/utils/paths.jl")

#auxiliary file directory definitions
const AUTOMATED_SIMULS_DIR_NAMES = readdir(AUTOMATED_SIMULS_DIR)

#= file manipulations =#
touch("median_magn_automated.txt") #creation of file containing all mean global magnetization
file_to_write_over = "median_magn_automated.txt"

#= plotting each temperature group and appending it over the file median_magn_automated.txt =#
GraphTrazes.graph_and_write_over_file!(AUTOMATED_SIMULS_DIR_NAMES,AUTOMATED_SIMULS_DIR,file_to_write_over)

# ploting mean magnetization vs temperature
GraphTrazes.plot_median_magn("median_magn_automated.txt", "median_magn_automated.pdf")