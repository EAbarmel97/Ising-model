module IsingModel
export do_simulations

include("Ising.jl")
using .Ising: isingModel 

include("IsingMethods.jl")
using .IsingMethods: display, reset_stats, compute_energy_cell, update_energy, update_magnetization,update_ising_model ,randomize, set_magnetization
using .IsingMethods: get_cell_coords, get_cell_id, do_generation, choose_flip_strategy, set_flip_strategy_and_transition_dynamics

include("utils/utilities.jl")
using .utilities: create_simulations_dir_if_not_exist

function do_model(init_magn::Float64, temp::Float64, n_grid::Float64, write_evol_array::Bool=false, is_automated::Bool=true)
    if is_automated
        utilities.create_automated_simulations_dir_if_not_exists()
    else
        utilities.create_simulations_dir_if_not_exists()
    end
    
    ising_model = IsingMethods.isingModel(temp, n_grid) #ising model struct instantiation
    
    IsingMethods.set_flip_strategy_and_transition_dynamics(ising_model,is_automated)
    
    simulations_dir = utilities.create_simulation_sub_dir(temp,is_automated)
    
    fourier_dir = utilities.create_fourier_dir(simulations_dir)
    magnetization_dir = utilities.create_magnetization_dir(simulations_dir)
    
    for run in 1:ARGS[3]
        IsingMethods.update_ising_model(ising_model,init_magn)
        
        magnetization_file_path = utilities.create_magnetization_time_series_file(magnetization_dir,"global_magnetization_r$(run).txt")
        IsingMethods.write_ising_model_prop_initial_state_over_file(ising_model,magnetization_file_path,:global_magnetization)
             
        for generation in 1:ARGS[6]
            IsingMethods.do_generation_and_write_ising_model_prop_over_file(ising_model, magnetization_file_path,:global_magnetization,generation)
            
            if write_evol_array
                IsingMethods.write_ising_model_sprin_grid(ising_model,generic_magnetization_file_name, generation)
            end     
        end
    end
end

#= do_model function wrapper to make simulations of the ising model at different temperatures =#
function do_simulations(ARGS::Array{Union{Int64,Float64},1},is_automated::Bool=false)
    
    if ARGS[2] < ARGS[1]
        throw(Exceptions.IlegalChoiceException("Ilegal  choice. Tf < Ti"))   
    end

    #= number of different temperatures equaly spaced by given incrments, contained in the interval [Ti, Tf] =#
    num_temps = ceil(Int,(ARGS[1] - ARGS[1])/ARGS[4])
    for i in 0:num_temps
        #= random initial temperature on the interval [-1 ,1] =#
        rand_magn = rand()*2 - 1 
            
        #= temperature increments in arithmetic porgression  =#
        temp = ARGS[1] + (i/num_temps)*(ARGS[2]- ARGS[1])
            
        do_model(rand_magn, temp, ARGS[5]) 
    end   
end
end #end of module