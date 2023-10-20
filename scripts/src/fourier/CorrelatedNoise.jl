module CorrelatedNoise
export correlated_noise_generator,imexp

function imexp(x::Array{Float64,1})::Array{ComplexF64}
    unitary_norm_complex_arr = ComplexF64[]
    for i in eachindex(x)
        unitary_norm_complex= cos(x[i]) + im*sin(x[i])
        push!(unitary_norm_complex_arr,unitary_norm_complex)
    end
    return unitary_norm_complex_arr
end

function correlated_noise_generator(N::Int64,beta0::Float64,beta1::Float64)::Array{ComplexF64,1}
    index_arr = collect(1:N)
    log_psd = map(index_arr) do u
        return beta0 + beta1*u
    end

    psd = exp.(log_psd)
    norm_array = sqrt.(psd)

    phase_arr = 2*pi*rand(N)

    return norm_array .* CorrelatedNoise.imexp(phase_arr)
end

end #end of module