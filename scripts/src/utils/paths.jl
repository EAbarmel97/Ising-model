"""
path config definitions
"""

const REPO_DIR = abspath(joinpath(@__DIR__,"..", "..", ".."))

#= simulations =#
const SIMULS_DIR = joinpath(REPO_DIR,"all_simulations")
const AUTOMATED_SIMULS_DIR = joinpath(SIMULS_DIR, "automated")

#= graphs dir and sub_dirs =#
const GRAPHS_DIR = joinpath(REPO_DIR,"graphs")
const GRAPHS_DIR_SIMULS = joinpath(GRAPHS_DIR,"simuls")
const PSD_GRAPHS = joinpath(GRAPHS_DIR,"psd")
const PSD_GRAPHS_SIMULS = joinpath(PSD_GRAPHS,"simuls")
const AUTOMATED_GRAPHS_DIR = joinpath(GRAPHS_DIR, "automated")
const AUTOMATED_GRAPHS_DIR_SIMULS = joinpath(AUTOMATED_GRAPHS_DIR,"simuls")
const AUTOMATED_PSD_GRAPHS = joinpath(AUTOMATED_GRAPHS_DIR,"psd")
const AUTOMATED_PSD_GRAPHS_SIMULS = joinpath(AUTOMATED_PSD_GRAPHS,"simuls")

#= fourier analysis =#
const ALL_GLOBAL_MAGN_DIRS = joinpath.(AUTOMATED_SIMULS_DIR, readdir(AUTOMATED_SIMULS_DIR), "magnetization")
const ALL_AUTOMATED_RFFTS = joinpath.(AUTOMATED_SIMULS_DIR, readdir(AUTOMATED_SIMULS_DIR), "fourier")

#= correlated noise =#
const CORRELATED_NOISE_DIR = joinpath(AUTOMATED_GRAPHS_DIR,"correlated_noise")

#= eigen spectrum =#
const AUTOMATED_EIGEN_SEPCTRUM_GRAPHS_DIR = joinpath(AUTOMATED_GRAPHS_DIR, "svd")