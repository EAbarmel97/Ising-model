include("../src/GraphTrazes.jl")
using .GraphTrazes: save_traze

include("../src/fourier/FourierAnalysis.jl")
using .FourierAnalysis: plot_psd

include("../src/fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator

include("../src/utils/paths.jl")

const path_to_save = joindir(AUTOMATED_PSD_GRAPHS,)

time_series1 = CorrelatedNoise.correlated_noise_generator(1000,1.0,0.0)
GraphTrazes.save_traze(path_to_save,time_series1)
FourierAnalysis.plot_psd(time_series1,"")
#println(time_series1)

time_series2 = CorrelatedNoise.correlated_noise_generator(1000,1.0,1.0)
GraphTrazes.save_traze(path_to_save,time_series2)
FourierAnalysis.plot_psd(time_series2,"")
#println(time_series2)

time_series3 = CorrelatedNoise.correlated_noise_generator(1000,1.0,2.0)
GraphTrazes.save_traze(path_to_save,time_series3)
FourierAnalysis.plot_psd(time_series2,"")
#println(time_series3)