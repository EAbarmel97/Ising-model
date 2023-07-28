#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze, graph_and_write_over_file!

#auxiliary file directory definitions
const CURR_DIR = pwd()
const AUTOMATED_SIMULS = CURR_DIR * "/all_simulations/automated/"
const AUTOMATED_SIMULS_DIR_NAMES= readdir(AUTOMATED_SIMULS)

#= regex'es =# 
rgx1 = r"T_1_\d{1,2}"
rgx2 = r"T_2_\d{1,2}"
rgx3 = r"T_3_\d{1,2}"

#= file manipulations =#
touch("mean_magn_automated.txt") #creation of file containing all mean global magnetization
file_to_write_over = CURR_DIR * "/mean_magn_automated.txt"

#writing headers over file
file = open("mean_magn_automated.txt","w+")
write(file,"temp,mean_magn\n")
close(file)

#= plotting each temperature group and appending it over the file mean_magn_automated.txt =#
graphTrazes.graph_and_write_over_file!(AUTOMATED_SIMULS_DIR_NAMES,AUTOMATED_SIMULS,file_to_write_over,rgx1)
graphTrazes.graph_and_write_over_file!(AUTOMATED_SIMULS_DIR_NAMES,AUTOMATED_SIMULS,file_to_write_over,rgx2)
graphTrazes.graph_and_write_over_file!(AUTOMATED_SIMULS_DIR_NAMES,AUTOMATED_SIMULS,file_to_write_over,rgx3)
