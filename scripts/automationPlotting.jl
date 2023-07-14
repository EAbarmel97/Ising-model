#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze, graph_and_write_over_file

#= file manipulations =#
cd("..")
touch("mean_magn_automated.txt") #creation of file containing all mean global magnetization
curr_dir = pwd()
file_to_write_over = curr_dir * "/mean_magn_automated.txt"

#writing headers over file
file = open("mean_magn.txt","w+")
write(file,"temp,mean_val")
close(file)

#auxiliary file directory definitions
curr_dir = pwd()
automated_simuls = curr_dir * "/all_simulations/automated/"
automated_simuls_dir_names = readdir(automated_simuls)

#= regex'es =# 
rgx1 = r"T_1.\d{2}"
rgx2 = r"T_2.\d{2}"
rgx3 = r"T_3.\d{2}"

#= plotting each temperature group and appending it over the file mean_magn_automated.txt =#
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx1)
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx2)
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx3)
