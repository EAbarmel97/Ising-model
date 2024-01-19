include("../src/svd/SVD.jl")
using .SVD: create_ts_matrix, centralize_matrix, plot_eigen_spectrum

include("../src/utils/paths.jl")

ts_matrix = SVD.create_ts_matrix(2.0,10000,10000)
M = SVD.centralize_matrix(ts_matrix)
SVD.plot_eigen_spectrum(AUTOMATED_EIGEN_SEPCTRUM_GRAPHS_DIR,M)