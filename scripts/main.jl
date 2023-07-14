include("../scripts/src/isingMethods.jl")
using .isingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS
using .isingMethods: display, reset_stats, compute_energy_cell, update_energy, update_magnetization, randomize, set_magnetization
using .isingMethods: get_cell_coords, get_cell_id, do_generation, choose_flip_strategy

include("../scripts/src/utilities.jl")
using .utilities: parse_int_float64, get_array_from_txt

println()

println("Number of different temperature to simulate")
const NUM_TEMPS = utilities.parse_int_float64(Int, readline())

println("Number of runs")
const NUM_RUNS = utilities.parse_int_float64(Int, readline())

println()

println("Number of generations")
const NUM_GENERATIONS = utilities.parse_int_float64(Int, readline())

println()

println("Grid size")
const N_GRID = utilities.parse_int_float64(Int, readline())

function do_model(INIT_MAGN, TEMP, N_GRID)
   cd("../all_simulations") #going up in the working directory
   curr_dir = pwd()
   ising_model = isingMethods.isingModel(TEMP, N_GRID) #ising model struct instantiation

   #= User chooses flip strategy and transition dynamics =#
   isingMethods.choose_flip_strategy(ising_model)
   isingMethods.choose_trans_dynamics(ising_model)

   if TEMP == isingMethods.CRITICAL_TEMP
      str_temp = replace("$(TEMP)", "2.269185314213020" => "Tc") 
   else
      str_temp = replace("$(TEMP)", "." => "_") #stringified temperature with "." replaced by "_"
   end       

   #= aux_dir = "../scripts/simulations_T_" * str_temp #folder containing simulations al temp str_temp  =#
   aux_dir = "simulations_T_" * str_temp #folder containing simulations al temp str_temp 

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
      write(generic_spin_grid_file, "$(stringified_grid_spin)\n")
      close(generic_spin_grid_file)

      for generation in 1:NUM_GENERATIONS
         isingMethods.do_generation(ising_model)
         setfield!(ising_model, :cur_gen, generation)

         generic_magnetization_file = open(generic_magnetization_file_name, "a+")
         write(generic_magnetization_file, "$(ising_model.global_magnetization)\n") #global magnetization observation at generation i 
         close(generic_magnetization_file)

         generic_spin_grid_file = open(generic_spin_grid_file_name, "a+")
         stringified_grid_spin = isingMethods.display(ising_model, ising_model.cur_gen)
         write(generic_spin_grid_file, "$(stringified_grid_spin)") #spin grid observation at generation i 
         close(generic_spin_grid_file)
      end
   end
end

function main()
   for i in 1:NUM_TEMPS
      println("Initial magnetization")
      INIT_MAGN = utilities.parse_int_float64(Float64, readline())

      println()

      println("Initial temperature")
      TEMP = utilities.parse_int_float64(Float64, readline())

      println()

      do_model(INIT_MAGN, TEMP, N_GRID)

      sleep(3)
end

main()