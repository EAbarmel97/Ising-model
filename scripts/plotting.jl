#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze, graph_and_write_over_file

#= file manipulations =#
cd("..")
touch("mean_magn.txt")
curr_dir = pwd()
file_to_write_over = curr_dir * "/mean_magn.txt"

#writing headers over file
file = open("mean_magn.txt","w+")
write(file,"temp,mean_val")
close(file)

#auxiliary file directory definitions
curr_dir = pwd()
all_simuls = curr_dir * "/all_simulations/"
simuls_dir_names = readdir(all_simuls)

#= regex'es =# 
rgx1 = r"T_1.\d{1}"
rgx2 = r"T_2.\d{1}"
rgx3 = r"T_Tc"
rgx4 = r"T_3.\d{1}"

#= plotting each temperature group and appending it over the file mean_magn.txt =#
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx1)
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx2)
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx3)
graphTrazes.graph_and_write_over_file(automated_simuls_dir_names,file_to_write_over,rgx4)
