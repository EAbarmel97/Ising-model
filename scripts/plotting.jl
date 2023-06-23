#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze

cd("..") #change current directory

#auxiliary file directory definitions
curr_dir = pwd()

file_dir_template1 = curr_dir * "/all_simulations/simulations_T_1"
file_dir_template2 = curr_dir * "/all_simulations/simulations_T_2"
file_dir_Tc = curr_dir * "/all_simulations/simulations_T_2_269185314213020/magnetizaton/global_magnetization_r1.txt"
file_dir_template3 = curr_dir * "/all_simulations/simulations_T_3"

for i in 0:9 #for temperatures in the interval 1.0-1.9 
    aux_dir = file_dir_template1 * "_$(i)" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = curr_dir * "/graphs/magn_ts_T_1_$(i).pdf"
    if isfile(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

for i in 0:9 #for temperatures in the interval 2.0-2.9
    aux_dir = file_dir_template2 * "_$i" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = curr_dir * "/graphs/magn_ts_T_1_$(i).pdf"
    if isfile(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

#at critical temperature
if isfile(file_dir_Tc)
    graph_dir = curr_dir * "/graphs/magn_ts_Tc.pdf"
    graphTrazes.save_traze(graph_dir, file_dir_Tc)
end

for i in 0:5 #for temperatures in the interval 3.0-3.5
    aux_dir = file_dir_template3 * "_$(i)" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = curr_dir * "/graphs/magn_ts_T_1_$(i).pdf"
    if isfile(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

