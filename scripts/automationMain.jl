include("../scripts/src/isingMethods.jl")
using .isingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS
using .isingMethods: display, reset_stats, compute_energy_cell, update_energy, update_magnetization, randomize, set_magnetization
using .isingMethods: get_cell_coords, get_cell_id, do_generation, choose_flip_strategy

include("../scripts/src/utilities.jl")
using .utilities: parse_int_float64, get_array_from_txt

using Random
Random.seed!(1234)

println("Provide initial and final temperatures. Ti < Tf")

println()

println("Initial temperature:")
const INIT_TEMP = utilities.parse_int_float64(Float64, readline())

println()

println("Final temperature:")
const FINAL_TEMP = utilities.parse_int_float64(Float64, readline())

println()

println("Number of runs:")
const NUM_RUNS = utilities.parse_int_float64(Int, readline())

println()

println("Number of temperatures")
const NUM_TEMPS = utilities.parse_int_float64(Int, readline())

println()

println("Grid size")
const N_GRID = utilities.parse_int_float64(Int, readline())

println()

println("Number of generations:")
const NUM_GENERATIONS = utilities.parse_int_float64(Int, readline())

cd("all_simulations/automated") #going up in the working directory

function do_model(INIT_MAGN, TEMP, N_GRID)
   curr_dir = pwd()
   println("saving simulations under dir: $curr_dir")
   ising_model = isingMethods.isingModel(TEMP, N_GRID) #ising model struct instantiation

   ising_model.flip_strategy = isingMethods.RANDOM_STRATEGY
   ising_model.trans_dynamics = isingMethods.METROPOLIS_DYNAMICS
   
   ROUNDED_TEMP = round(TEMP, digits=2)
   str_temp = replace("$(ROUNDED_TEMP)", "." => "_") #stringified temperature with "." replaced by "_"
   aux_dir =  curr_dir * "/simulations_T_" * str_temp #folder containing simulations al temp str_temp 

   mkpath(aux_dir) #creates simulation folder  

   #= Global magnetization time series realization will be saved on subdirectories over folder simultations=#
   global_magnetization_aux_dir = aux_dir * "/magnetization"
   mkpath(global_magnetization_aux_dir)

   #= Subdirectory containg a .txt file with the unicode representation of how the spin grid evolves with each generation at each run =#
   grid_evolution_aux_dir = aux_dir * "/grid_evolution"
   mkpath(grid_evolution_aux_dir)

   for run in 1:NUM_RUNS
      isingMethods.reset_stats(ising_model)
      isingMethods.set_magnetization(INIT_MAGN, ising_model) #populates the spin grid with a given initial magnetization 
      isingMethods.update_magnetization(ising_model) #updates global magnetization 
      isingMethods.update_energy(ising_model) #updates global energy

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
      stringified_grid_spin = isingMethods.display(ising_model, ising_model.cur_gen)
      write(generic_spin_grid_file, "$(stringified_grid_spin)")
      close(generic_spin_grid_file)

      for generation in 1:NUM_GENERATIONS
         isingMethods.do_generation(ising_model)
         setfield!(ising_model, :cur_gen, generation)

         generic_magnetization_file = open(generic_magnetization_file_name, "a+")
         write(generic_magnetization_file, "$(ising_model.global_magnetization)\n") #global magnetization observation at generation i 
         close(generic_magnetization_file)

         generic_spin_grid_file = open(generic_spin_grid_file_name, "a+")
         stringified_grid_spin = isingMethods.display(ising_model, ising_model.cur_gen)
         write(generic_spin_grid_file, "$(stringified_grid_spin)\n") #spin grid observation at generation i 
         close(generic_spin_grid_file)
      end
   end
end

function main()
   for i in 0:NUM_TEMPS
      #= random initial temperature on the interval [-1 ,1] =#
      RAND_MAGN = rand()*2 - 1 
      
      #= temperature increments in arithmetic porgression  =#
      temp = INIT_TEMP + (i/NUM_TEMPS)*(FINAL_TEMP - INIT_TEMP)

      do_model(RAND_MAGN, temp, N_GRID) 

      sleep(3)
   end
end

main()