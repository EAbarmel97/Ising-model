module CorrelatedNoise
using FFTW

export correlated_noise_generator,imexp

"""
    imexp(x::Array{Float64,1})::Array{ComplexF64,1}
    
Returns an array of unitary norm complex numbers. Obtained from Euler's polar form
"""
function imexp(x::Array{Float64,1})::Array{ComplexF64,1}
    unitary_norm_complex_arr = ComplexF64[]
    for i in eachindex(x)
        unitary_norm_complex = cos(x[i]) + sin(x[i])*im
        push!(unitary_norm_complex_arr,unitary_norm_complex)
    end
    return unitary_norm_complex_arr
end

"""
    correlated_noise_generator(N::Int64,beta0::Float64,beta1::Float64)::Array{Float64,1} 

Outputs a time series obtained by inverse-fourier transforming correlated noise given exponenet and intercept from a linear fit 
"""
function correlated_noise_generator(N::Int64,beta0::Float64,beta1::Float64)::Array{Float64,1}
    index_arr = collect(1:N)
    log10_psd = map(index_arr) do u
        return beta0 - beta1*u
    end

    psd = exp10.(log10_psd)
    norm_array = sqrt.(psd) 
    phase_arr = 2*pi*rand(N) 
    z_array = norm_array .* imexp(phase_arr) #correlated noise

    return irfft(z_array,2*(length(z_array) - 1))
end
end #end of module