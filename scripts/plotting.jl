#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.2=#

include("../scripts/.julia/src/graphTrazes.jl")
using .graphTrazes: save_traze

file_dir_template1 = "../scripts/simulations_T_1"
file_dir_template2 = "../scripts/simulations_T_2"
file_dir_template3 = "../scripts/simulations_T_3"

for i in 1:9 #for temperatures in the interval 1.0-1.9 
    aux_dir *= file_dir_template1 * "_$i" * "/magnetization/" * "magn_ts_$(i).pdf"
    graphTrazes.save_traze(aux_dir, "") 
end

for i in 1:9 #for temperatures in teh interval 2.1
    aux_dir *= file_dir_template2 * "_$i" * "/magnetization/" * "magn_ts_$(i).pdf"
    graphTrazes.save_traze(aux_dir,"") 
end

for i in 1:3 #for temperatures in the interval 
    aux_dir *= file_dir_template3 * "_$i" * "/magnetization/" * "magn_ts_$(i).pdf"
    graphTrazes.save_traze(aux_dir,"") 
end