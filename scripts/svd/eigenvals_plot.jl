include("../src/svd/SVD.jl")
using .SVD: create_ts_matrix, centralize_matrix, plot_eigen_spectrum, write_beta_beta_fit, plot_beta_beta_fit

include("../src/utils/paths.jl")

ts_matrix = SVD.create_ts_matrix(2.0;number_of_observations=1000,number_of_realizations=10000)
M = SVD.centralize_matrix(ts_matrix)
SVD.plot_eigen_spectrum(AUTOMATED_EIGEN_SEPCTRUM_GRAPHS_DIR,M,2.0)

write_beta_beta_fit(0.0,3.0,2;number_of_observations=100,number_of_realizations=100)
plot_beta_beta_fit("beta_beta_fit.txt")