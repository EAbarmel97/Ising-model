module Ising
export isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS 
#The critical temperature is Tc = 2/ln(1+sqrt(2)) in units of J/k, with k beeing the Boltzman constant
const CRITICAL_TEMP = 2.26918531421302

const RANDOM_STRATEGY = 0
const SHUFFLE_STRATEGY = 1
const SEQUENTIAL_STRATEGY = 2

const METROPOLIS_DYNAMICS = 0
const GLAUBER_DYNAMICS = 1

mutable struct isingModel
    #Temperature 
    TEMP :: Float64

    #Size of grid 
    NGRID :: Int

    #Size of grid's side
    NCELLS :: Int

    #2D array containg the spin states: N²
    grid :: Array{Int,2}

    #Flip order 
    flip_order :: Array{Int,1}

    #Fliping strategies 
    flip_strategy :: Int
    
    #Transition dynamics 
    trans_dynamics :: Int 
    
    #Current generation (will never reset)
    cur_gen :: Int 
    
    #Global statistics
    global_energy :: Float64
    global_mean :: Float64
    global_variance :: Float64
    global_magnetization :: Float64

    #model constructor with two params 
    function isingModel(p_TEMP,p_NGRID)
        TEMP = p_TEMP
        NGRID = p_NGRID
        NCELLS = NGRID*NGRID
        flip_strategy = RANDOM_STRATEGY
        flip_order = Array{Int,1}(undef,NCELLS)
        trans_dynamics = METROPOLIS_DYNAMICS
        grid = Array{Int,2}(undef,NGRID,NGRID)
        cur_gen = 0 

        global_energy = 0.0
        global_magnetization = 0.0
        global_mean = 0.0
        global_variance = 0.0

        return new(TEMP,NGRID,NCELLS,grid,flip_order,flip_strategy,trans_dynamics,cur_gen, 
        global_energy, global_magnetization,global_mean, global_variance)
    end
    
end #end of struct 
end
