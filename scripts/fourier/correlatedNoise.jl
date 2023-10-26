include("../src/GraphTrazes.jl")
using .GraphTrazes: save_traze

include("../src/fourier/FourierAnalysis.jl")
using .FourierAnalysis: plot_psd

include("../src/fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator

include("../src/utils/utilities.jl")
using .utilities: determines_noise_or_movement

include("../src/utils/paths.jl")

function create_correlated_noise_dir()
    if !isdir(CORRELATED_NOISE_DIR)
        mkpath(CORRELATED_NOISE_DIR)
    end
end

function correlated_noise_graph_file_path(dir::String,beta0::Float64,beta1::Float64)::String
    return utilities.psd_graph_file_path(dir,beta0,beta1)  
end

create_correlated_noise_dir()

time_series1 = CorrelatedNoise.correlated_noise_generator(1001,1.0,0.0)
path_to_save = correlated_noise_graph_file_path(CORRELATED_NOISE_DIR,1.0,0.0)
if !isfile(path_to_save) 
    GraphTrazes.save_traze(path_to_save,time_series1,1.0,0.0)    
    FourierAnalysis.plot_psd(time_series1,CORRELATED_NOISE_DIR,1.0,0.0)
end

time_series2 = CorrelatedNoise.correlated_noise_generator(1001,1.0,1.0)
path_to_save = correlated_noise_graph_file_path(CORRELATED_NOISE_DIR,1.0,1.0)
if !isfile(path_to_save)
    GraphTrazes.save_traze(path_to_save,time_series2,1.0,1.0)
    FourierAnalysis.plot_psd(time_series2,CORRELATED_NOISE_DIR,1.0,1.0)    
end

time_series3 = CorrelatedNoise.correlated_noise_generator(1001,1.0,2.0)
path_to_save = correlated_noise_graph_file_path(CORRELATED_NOISE_DIR,1.0,2.0)
if !isfile(path_to_save)
    GraphTrazes.save_traze(path_to_save,time_series3,1.0,2.0)
    FourierAnalysis.plot_psd(time_series2,CORRELATED_NOISE_DIR,1.0,2.0)
end
