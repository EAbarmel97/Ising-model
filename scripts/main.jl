include("../scripts/.julia/src/isingMethods.jl")
using .isingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS,GLAUBER_DYNAMICS
using .isingMethods: display,reset_stats,compute_energy_cell,update_energy,update_magnetization,randomize,set_magnetization
using .isingMethods: get_cell_coords, get_cell_id, do_generation,choose_flip_strategy

include("../scripts/.julia/src/utilities.jl")
using .utilities: parse_int_float64

println("provide an initial magnetization")
const INIT_MAGN = utilities.parse_int_float64(Float64,readline())

println()

println("provide the number of runs")
const NUM_RUNS = utilities.parse_int_float64(Int,readline())

println()

println("provide the number of generations")
const NUM_GENERATIONS = utilities.parse_int_float64(Int,readline()) 

println()

println("provide an initial temperature")
const TEMP = utilities.parse_int_float64(Float64,readline())

println()

println("provide a grid size")
const N_GRID = utilities.parse_int_float64(Int,readline())

println()

function main()

   ising_model = isingMethods.isingModel(TEMP,N_GRID) #ising model struct instantiation 
   
   #= User chooses flip strategy and transition dynamics =#
   isingMethods.choose_flip_strategy(ising_model)
   isingMethods.choose_trans_dynamics(ising_model)

   cd("scripts") #cd to scripts subdirectory
 
   #= Global energy and magnetization time series realization will be saved on subdirectories over folder simultations=#
   mkpath("simulations/energy")
   mkpath("simulations/magnetization")

   #= Subdirectories containg a .txt file with the unicode representation of how the spin grid evolves with each generation at each run =#
   mkpath("simulations/grid_evolution")
   
   file_names = ["global_energy","global_magnetization"] 

   for run in 1:NUM_RUNS      
      isingMethods.reset_stats(ising_model) 
      isingMethods.set_magnetization(INIT_MAGN,ising_model) #populates the spin grid with a given initial magnetization 
      isingMethods.update_magnetization(ising_model) #updates global magnetization 
      isingMethods.update_energy(ising_model) #updates global energy
      
      #= Cration of generic .txt files containing energy and magnetization time series =#
      generic_energy_file_name = "simulations/energy/$(file_names[1])_r$(run)"*".txt"
      touch("$(generic_energy_file_name)")
      
      generic_magnetization_file_name = "simulations/magnetization/$(file_names[2])_r$(run)"*".txt"
      touch("$(generic_magnetization_file_name)") 
      
      #= Cration of generic .txt files containing snapshots of the spin grid evolution at each generation =#
      generic_spin_grid_file_name = "simulations/grid_evolution/grid_evolution_r$(run)"*".txt"
      touch("$(generic_spin_grid_file_name)")

      #= Initial observations of the global energy and and magnetizaton are saved to their respective .txt files=#
      generic_energy_file = open(generic_energy_file_name,"a+")
      write(generic_energy_file_name, "$(ising_model.global_energy)\n")
      close(generic_energy_file)
      
      generic_magnetization_file = open(generic_magnetization_file_name,"a+")
      write(generic_energy_file_name, "$(ising_model.global_magnetization)\n")
      close(generic_energy_file)

      #= Initial spin grid state =#
      generic_spin_grid_file = open(generic_spin_grid_file_name,"a+")
      stringified_grid_spin = isingMethods.display(ising_model,ising_model.cur_gen)
      write(generic_spin_grid_file,"$(stringified_grid_spin)\n")
      close(generic_spin_grid_file)
      
      for generation in 1:NUM_GENERATIONS 
         isingMethods.do_generation(ising_model) 
         setfield!(ising_model,:cur_gen, generation)

         #= Writing observation i to the global energy  and magnetization time series=# 
         generic_energy_file = open(generic_energy_file_name,"a+")
         write(generic_energy_file_name, "$(ising_model.global_energy)\n") #global energy observation at generation i 
         close(generic_energy_file)

         generic_magnetization_file = open(generic_magnetization_file_name,"a+")
         write(generic_magnetization_file,"$(ising_model.global_magnetization)\n") #global magnetization at generation i 
         close(generic_magnetization_file)
         
         generic_spin_grid_file = open(generic_spin_grid_file_name,"a+")
         stringified_grid_spin = isingMethods.display(ising_model,ising_model.cur_gen)
         write(generic_spin_grid_file,"$(stringified_grid_spin)") #spin grid observation at generation i 
         close(generic_spin_grid_file)
      end      
   end
   
   #= TO DO: implement graphs of global energy and global magnetization time series =#

end

main()





