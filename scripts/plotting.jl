#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze, graph_and_write_over_file!

#auxiliary file directory definitions
const CURR_DIR = pwd()
const ALL_SIMULS = CURR_DIR * "/all_simulations/"
const SIMULS_DIR_NAMES = readdir(ALL_SIMULS)

#= regex'es =# 
rgx1 = r"T_1_\d{1}"
rgx2 = r"T_2_\d{1}"
rgx3 = r"T_Tc"
rgx4 = r"T_3_\d{1}"

#= file manipulations =#
touch("mean_magn.txt")
file_to_write_over = CURR_DIR * "/mean_magn.txt"

#writing headers over file
file = open("mean_magn.txt","w+")
write(file,"temp,mean_val\n")
close(file)

#= plotting each temperature group and appending it over the file mean_magn.txt =#
#graphTrazes.graph_and_write_over_file!(SIMULS_DIR_NAMES,ALL_SIMULS,file_to_write_over,rgx1)
graphTrazes.graph_and_write_over_file!(SIMULS_DIR_NAMES,ALL_SIMULS,file_to_write_over,rgx2)
#graphTrazes.graph_and_write_over_file!(SIMULS_DIR_NAMES,ALL_SIMULS,file_to_write_over,rgx3)
#graphTrazes.graph_and_write_over_file!(SIMULS_DIR_NAMES,ALL_SIMULS,file_to_write_over,rgx4)
