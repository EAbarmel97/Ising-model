include("../src/svd/SVD.jl")
using .SVD: create_ts_matrix, centralize_matrix, plot_eigen_spectrum, write_beta_beta_fit, read_beta_beta_fit, plot_beta_beta_fit, average_eigen_spectrum

include("../src/utils/paths.jl")

function write_beta_beta_fit_cli(ARGS)
    arg1 = Base.parse(Int64,ARGS[1])
    arg2 = Base.parse(Int64,ARGS[2])
    write_beta_beta_fit(0.0,3.0,2;number_of_observations=arg1,number_of_realizations=arg2)
end

function plot_average_eigenspec_cli(ARGS)
    arg1 = Base.parse(Int64,ARGS[1])
    arg2 = Base.parse(Int64,ARGS[2])
    println(average_eigen_spectrum(2.0,number_of_observations=arg1,number_of_realizations=arg2, num_samples=10))
end

write_beta_beta_fit_cli(ARGS)
plot_average_eigenspec_cli(ARGS)
plot_beta_beta_fit(read_beta_beta_fit, "beta_beta_fit.txt")