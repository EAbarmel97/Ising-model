using Random
Random.seed!(1234)

include("../scripts/src/IsingMethods.jl")
using .IsingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS
using .IsingMethods: display, reset_stats, compute_energy_cell, update_energy, update_magnetization, randomize, set_magnetization
using .IsingMethods: get_cell_coords, get_cell_id, do_generation, choose_flip_strategy

include("../scripts/src/utils/utilities.jl")
using .utilities: parse_int_float64, get_array_from_txt, parse_int_float64, use_temperature_array, TEMPERATURE_INTERVALS

const SIMULS_DIR = "all_simulations"

if !isdir(SIMULS_DIR)
   mkpath(SIMULS_DIR)
end

println()

println("Number of different temperature to simulate\n")
const NUM_TEMPS = utilities.parse_int_float64(Int, readline())

println()

println("Number of runs\n")
const NUM_RUNS = utilities.parse_int_float64(Int, readline())

println()

println("Number of generations\n")
const NUM_GENERATIONS = utilities.parse_int_float64(Int, readline())

println()

println("Grid size\n")
const N_GRID = utilities.parse_int_float64(Int, readline())

cd("all_simulations") #going up in the working directory

const CURR_DIR = pwd()

println()

println("saving simulations under dir: $CURR_DIR\n") 

function do_model(INIT_MAGN, TEMP, N_GRID)
   ising_model = IsingMethods.isingModel(TEMP, N_GRID) #ising model struct instantiation

   #= User chooses flip strategy and transition dynamics =#
   IsingMethods.choose_flip_strategy(ising_model)
   IsingMethods.choose_trans_dynamics(ising_model)

   if TEMP == IsingMethods.CRITICAL_TEMP
      str_temp = replace("$(TEMP)", "2.269185314213020" => "Tc") 
   else
      str_temp = replace("$(TEMP)", "." => "_") #stringified temperature with "." replaced by "_"
   end       

   #= aux_dir = "../scripts/simulations_T_" * str_temp #folder containing simulations al temp str_temp  =#
   aux_dir = CURR_DIR * "/simulations_T_" * str_temp #folder containing simulations al temp str_temp 
   FOURIER_DIR = joinpath(aux_dir,"fourier")
   
   mkpath(aux_dir) #creates simulation folder
   mkpath(FOURIER_DIR)

   #= Global magnetization time series realization will be saved on subdirectories over folder simultations=#
   global_magnetization_aux_dir = joinpath(aux_dir, "magnetization") 
   mkpath(global_magnetization_aux_dir)

   #= Subdirectory containg a .txt file with the unicode representation of how the spin grid evolves with each generation at each run =#
   grid_evolution_aux_dir = aux_dir * "/grid_evolution"
   mkpath(grid_evolution_aux_dir)

   for run in 1:NUM_RUNS
      IsingMethods.reset_stats(ising_model)
      IsingMethods.set_magnetization(INIT_MAGN, ising_model) #populates the spin grid with a given initial magnetization 
      IsingMethods.update_magnetization(ising_model) #updates global magnetization 
      IsingMethods.update_energy(ising_model) #updates global energy

      #= Creation of generic .txt files containing global magnetization time series =#
      generic_magnetization_file_name = global_magnetization_aux_dir * "/global_magnetization_r$(run)" * ".txt"
      touch("$generic_magnetization_file_name")

      #= Creation of generic .txt files containing snapshots of the spin grid evolution at each generation =#
      generic_spin_grid_file_name = grid_evolution_aux_dir * "/grid_evolution_r$(run)" * ".txt"
      touch("$generic_spin_grid_file_name")

      #= Initial observations of the global magnetizaton are saved to their respective .txt files=#
      generic_magnetization_file = open(generic_magnetization_file_name, "w+")
      write(generic_magnetization_file, "$(ising_model.global_magnetization)\n") #initial observation
      close(generic_magnetization_file)

      #= Initial spin grid state =#
      generic_spin_grid_file = open(generic_spin_grid_file_name, "w+")
      stringified_grid_spin = IsingMethods.display(ising_model, ising_model.cur_gen)
      write(generic_spin_grid_file, "$(stringified_grid_spin)")
      close(generic_spin_grid_file)

      for generation in 1:NUM_GENERATIONS
         IsingMethods.do_generation(ising_model)
         setfield!(ising_model, :cur_gen, generation)

         generic_magnetization_file = open(generic_magnetization_file_name, "a+")
         write(generic_magnetization_file, "$(ising_model.global_magnetization)\n") #global magnetization observation at generation i 
         close(generic_magnetization_file)

         generic_spin_grid_file = open(generic_spin_grid_file_name, "a+")
         stringified_grid_spin = IsingMethods.display(ising_model, ising_model.cur_gen)
         write(generic_spin_grid_file, "$(stringified_grid_spin)\n") #spin grid observation at generation i 
         close(generic_spin_grid_file)

         if generation == NUM_GENERATIONS
            IsingMethods.do_generation(ising_model)
            setfield!(ising_model, :cur_gen, NUM_GENERATIONS)
      
            generic_magnetization_file = open(generic_magnetization_file_name, "a+")
            write(generic_magnetization_file, "$(ising_model.global_magnetization)") #global magnetization observation at the last generation
            close(generic_magnetization_file)
      
            generic_spin_grid_file = open(generic_spin_grid_file_name, "a+")
            stringified_grid_spin = IsingMethods.display(ising_model, ising_model.cur_gen)
            write(generic_spin_grid_file, "$(stringified_grid_spin)") #spin grid observation at the last generation
            close(generic_spin_grid_file)  
         end
      end     
   end
end

function do_simulations(arr :: Array{Float64,1})
   for i in eachindex(arr)
      #= random initial temperature on the interval [-1 ,1] =#
      rand_magn = rand()*2 - 1 
      
      #= temperature increments in arithmetic progression  =#
      temp = arr[i]

      do_model(rand_magn, temp, N_GRID) 
   end
end

#= 
   if user wants to simulate with the default array of temperatures containing temps in the intervals
   a) 0.0 to 1.0 with increments of 0.1
   b) 1.0 to 2.2 with increments of 0.1 
   c) 2.2 to 2.5 with incremets of 0.01 
   d) 2.5 to 3.5 with increments of 0.1 
=#

function main()
   println()
   if utilities.use_temperature_array()
      do_simulations(utilities.TEMPERATURE_INTERVALS)
   else
      for i in 1:NUM_TEMPS
         println()

         println("Initial magnetization")
         INIT_MAGN = utilities.parse_int_float64(Float64, readline())
   
         println()
   
         println("Initial temperature")
         TEMP = utilities.parse_int_float64(Float64, readline())
   
         println()
   
         do_model(INIT_MAGN, TEMP, N_GRID)
      end
   end     
end

main()