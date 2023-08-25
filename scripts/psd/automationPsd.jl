include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd, sampling_freq_arr

include("../src/utilities.jl")
using .utilities: get_array_from_txt, parse_int_float64

#= auxiliary constants  =#
const CURR_DIR = pwd()
const AUTOMATED_SIMULS_DIR = CURR_DIR * "/all_simulations/automated/"
const AUTOMATED_PSD_GRAPHS = "graphs/automated/psd/"
const ALL_AUTOMATED_SIMULS_DIRS = readdir(AUTOMATED_SIMULS_DIR)
const ALL_GLOBAL_MAGN_DIRS = string.(AUTOMATED_SIMULS_DIR,ALL_AUTOMATED_SIMULS_DIRS,"/magnetization/")

#number of runs equals the number of files in ../magnetization/ and is the same for each simulations_T_x_y_z/magnetization dir 
const NUM_RUNS = length(readdir(ALL_GLOBAL_MAGN_DIRS[1]))

#create dir if it doesn't already exist
if !isdir(AUTOMATED_PSD_GRAPHS)
    mkpath(AUTOMATED_PSD_GRAPHS)
end

for run in 1:NUM_RUNS
    for i in eachindex(ALL_AUTOMATED_SIMULS_DIRS)
        simul_dir_name = ALL_AUTOMATED_SIMULS_DIRS[i]
    
        #global magnetization with initial temperature T_x_y_z
        global_magn_ts_path = simul_dir_name * readdir(simul_dir_name)[1]
    
        #temperature is taken from simulations dir name 
        str_temp = replace(simul_dir_name,"/all_simulations/automated/simulations_T_" => "", "_" => ".")
        
        #dir where rfft will be saved
        fourier_dir = ALL_AUTOMATED_SIMULS_DIRS[i] * "/fourier/"

        #if strigified rfft file doesn't exist at dir ../fourier/ 
        if isempty(readdir(fourier_dir))
            #rfft is computed from .txt files containing the global magnetization time series
            rfft = fourierAnalysis.compute_rfft(global_magn_ts_path)

            temp = utilities.parse_int_float64(Float64,str_temp)

            #rfft is saved  as a .txt file under the dir ../automated/fourier/
            fourierAnalysis.write_rfft(rfft,fourier_dir,temp)
        end                                                             
    end    
end

   