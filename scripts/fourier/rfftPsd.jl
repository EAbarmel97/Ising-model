include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: compute_rfft, compute_psd, sampling_freq_arr, plot_psd, create_order_coef_dir_and_file

include("../src/utils/utilities.jl")
using .utilities: get_array_from_txt, parse_int_float64

include("../src/utils/paths.jl")

const ALL_AUTOMATED_SIMULS_DIRS = readdir(AUTOMATED_SIMULS_DIR)
#=
number of runs equals the number of files in ../magnetization/ and is the same for each simulations_T_x_y_z/magnetization dir
so the first simulation directory is choosen 
=#
const NUM_RUNS = length(readdir(ALL_GLOBAL_MAGN_DIRS[1]))

#create dir if it doesn't already exist
if !isdir(AUTOMATED_PSD_GRAPHS)
    mkpath(AUTOMATED_PSD_GRAPHS)
end

fourierAnalysis.create_order_coef_dir_and_file()

#writing under each simulations_T_x_y_z/fourier/ dir the rfft at each run and plotting the psd
for i in eachindex(ALL_AUTOMATED_SIMULS_DIRS)
    # array of strings has generic strings of the the type: simulations_T_x_y_z
    simul_dir_name = ALL_AUTOMATED_SIMULS_DIRS[i]
    simul_sub_dir = replace(ALL_AUTOMATED_SIMULS_DIRS[i],"simulations_" => "")

    for run in 1:NUM_RUNS
        #= 
        global magnetization with initial temperature T_x_y_z s is under the directory 
        ../automated/simulations_T_x_y_z/magnetization/global_magnetization_rW.txt 
        =#
        
        global_magn_ts_path = joinpath(ALL_GLOBAL_MAGN_DIRS[i],readdir(ALL_GLOBAL_MAGN_DIRS[i])[run])
        #temperature is taken from simulations dir name 
        str_temp = replace(simul_dir_name,"simulations_T_" => "", "_" => ".")
    
        #dir where rfft will be saved
        fourier_dir = joinpath(AUTOMATED_SIMULS_DIR, ALL_AUTOMATED_SIMULS_DIRS[i], "fourier/")
        rfft_path = joinpath(fourier_dir, "rfft_global_magnetization_$(str_temp)_r$run.txt")

        #if strigified rfft file doesn't exist at dir ../automated/simulations_T_x_y_z/fourier/
        if !isfile(rfft_path)
            #rfft is computed from .txt files containing the global magnetization time series
            rfft = fourierAnalysis.compute_rfft(global_magn_ts_path)
        
            temp = utilities.parse_int_float64(Float64,str_temp)
        
            #rfft is saved  as a .txt in ../automated/simulations_T_x_y_z/fourier/
            fourierAnalysis.write_rfft(rfft,fourier_dir,temp,run)        
        end
    end

    psd_plot_file_name = "psd_$(simul_sub_dir)_r_1_$(NUM_RUNS).pdf"
    psd_plot_file_abs_path = joinpath(AUTOMATED_PSD_GRAPHS,psd_plot_file_name)

    if !isfile(psd_plot_file_abs_path)
        #plotting the power density spectra
        fourierAnalysis.plot_psd(simul_dir_name,AUTOMATED_PSD_GRAPHS)      
    end

    #= TO DO: implement logic to obtain the order coefficient β =#
end