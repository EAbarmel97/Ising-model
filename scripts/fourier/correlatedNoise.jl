include("../src/fourier/CorrelatedNoise.jl")
using .CorrelatedNoise: correlated_noise_generator

time_series = CorrelatedNoise.correlated_noise_generator(1000,-2.3,0.7)

