include("../src/fourier/fourierAnalysis.jl")
using .fourierAnalysis: plot_psd

include("../src/fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator

time_series = CorrelatedNoise.correlated_noise_generator(1000,1.0,-1.0)

