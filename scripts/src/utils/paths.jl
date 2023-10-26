"""
path config definitions
"""
#= simulations =#
const SIMULS_DIR = "all_simulations"
const AUTOMATED_SIMULS_DIR = joinpath(SIMULS_DIR, "automated")

#= graphs dir and sub_dirs =#
const GRAPHS_DIR = "graphs"
const PSD_GRAPHS = joinpath(GRAPHS_DIR,"psd")
const AUTOMATED_GRAPHS_DIR = joinpath(GRAPHS_DIR, "automated")
const AUTOMATED_PSD_GRAPHS = joinpath(AUTOMATED_GRAPHS_DIR,"psd")

#= fourier analysis =#
const ALL_GLOBAL_MAGN_DIRS = joinpath.(AUTOMATED_SIMULS_DIR, readdir(AUTOMATED_SIMULS_DIR), "magnetization")
const ALL_AUTOMATED_RFFTS = joinpath.(AUTOMATED_SIMULS_DIR, readdir(AUTOMATED_SIMULS_DIR), "fourier")

#= correlated noise =#
const CORRELATED_NOISE = joinpath(AUTOMATED_GRAPHS_DIR,"correlated_noise")