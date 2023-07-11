#= Automation of global magnetization time series plot file saving at temperatures ranging from 1.0 to 3.5=#
include("../scripts/src/graphTrazes.jl")
using .graphTrazes: save_traze

include("../scripts/src/utilities")
using .utilities: mean_value

cd("..") #change current directory. Current directory goes up in the file hierarchy 

touch("mean_temps.txt") #file containing all mean temperatures is created

#= writing headers of the comma separted file =#
file = open("mean_temps.txt","w+")
write(file,"temp,mean_val")
close(file)

#auxiliary file directory definitions
curr_dir = pwd()
dir_names = readdir()

#regex'es 
rgx1 = r"T_1.\d{2}"
rgx2 = r"T_2.\d{2}"
rgx3 = r"T_3.\d{2}"

temps_1_xx = filter(file_name -> contains(file_name, rgx1), dir_names) #temperatures of the type 1.xx
temps_2_xx = filter(file_name -> contains(file_name, rgx2), dir_names) #temperatures of the type 2.xx
temps_3_xx = filter(file_name -> contains(file_name, rgx3), dir_names) #temperatures of the type 3.xx

general_aux_dir_template = curr_dir * "/all_simulations/automated/simulations_T"

for i in eachindex(temps_1_xx)
    aux_dir = general_aux_dir_template * "1_$(i)" * "/magnetization/global_magnetization_r1.txt"
    mean_val = utilities.mean_value(aux_dir) #mean value of the global magnetization time series
    str_mean_val = "$(mean_val)"

    #= appending mean value of the global magn time series =#
    str_to_append = str_mean_val
    file = open("mean_temps.txt", "a+")
    write(file, str_to str_to_append) 
    close(file)

    aux_graph_name = curr_dir * "/graphs/magn_ts_T_1_$(i).pdf"
    if isfile(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

for i in eachindex(temps_2_xx)
    aux_dir = file_dir_template2 * "_$i" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = curr_dir * "/graphs/magn_ts_T_1_$(i).pdf"
    if isfile(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

for i in eachindex(temps_3_xx)#temperatures in the interval 3.0-3.5
    aux_dir = file_dir_template3 * "_$(i)" * "/magnetization/global_magnetization_r1.txt"
    aux_graph_name = curr_dir * "/graphs/magn_ts_T_1_$(i).pdf"
    if isfile(aux_dir)
        graphTrazes.save_traze(aux_graph_name, aux_dir)
    end
end

