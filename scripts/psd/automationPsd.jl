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
        # the array of strings has generic strings of the simulations_T_x_y_z
        simul_dir_name = ALL_AUTOMATED_SIMULS_DIRS[i]
    
        #= global magnetization with initial temperature T_x_y_z s is under the directory
        ../automated/simulations_T_x_y_z/magnetization/global_magnetization_rW.txt =#
        global_magn_ts_path =  ALL_GLOBAL_MAGN_DIRS[i] * (readdir(ALL_GLOBAL_MAGN_DIRS[i])[run])
    
        #temperature is taken from simulations dir name 
        str_temp = replace(simul_dir_name,"simulations_T_" => "", "_" => ".")
        
        #dir where rfft will be saved
        fourier_dir =  AUTOMATED_SIMULS_DIR * ALL_AUTOMATED_SIMULS_DIRS[i] * "/fourier/"
        rfft_path = fourier_dir * fourier_dir * "rfft_global_magnetization_$str_temp$run.txt"

        #if strigified rfft file doesn't exist at dir ../automated/simulations_T_x_y_z/fourier/
        if !isfile(rfft_path)
            #rfft is computed from .txt files containing the global magnetization time series
            rfft = fourierAnalysis.compute_rfft(global_magn_ts_path)

            temp = utilities.parse_int_float64(Float64,str_temp)
            
            #rfft is saved  as a .txt in ../automated/simulations_T_x_y_z/fourier/
            fourierAnalysis.write_rfft(rfft,fourier_dir,temp,run)
        end                                                             
    end 
    
    #= 
    TO DO: change logic disign to be able to interactively ask to the user how to plot the power denstiy spectra 
    =#
    
end

   