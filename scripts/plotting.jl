#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#

include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze

cd("../scripts")#

#auxiliary file directory definitions
file_dir_template1 = "simulations_T_1"
file_dir_template2 = "simulations_T_2"
file_dir_Tc = "simulations_T_2_269185314213020/magnetizaton/global_magnetization_r1.txt"
file_dir_template3 = "simulations_T_3"


for i in 0:9 #for temperatures in the interval 1.0-1.9 
    aux_dir *= file_dir_template1 * "_$(i)" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = "graphs/" * "magn_ts_T_1_$(i).pdf"
    if isdir(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

for i in 0:9 #for temperatures in the interval 2.0-2.9
    aux_dir = file_dir_template2 * "_$i" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = "graphs/" * "magn_ts_T_1_$(i).pdf"
    if isdir(aux_dir)
        graphTrazes.save_traze("", aux_dir)
    end
end

#at critical temperature
if isdir(file_dir_Tc)
    graphTrazes.save_traze("graphs/magn_ts_Tc.pdf", file_dir_Tc)
end

for i in 0:5 #for temperatures in the interval 3.0-3.5
    aux_dir = file_dir_template3 * "_$(i)" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = "graphs/" * "magn_ts_T_1_$(i).pdf"
    if isdir(aux_dir)
        graphTrazes.save_traze("", aux_dir)
    end
end

graphTrazes.save_traze("graphs/magn_ts1.pdf", "scripts/simulations_T_1_0/magnetization/global_magnetization_r1.txt")