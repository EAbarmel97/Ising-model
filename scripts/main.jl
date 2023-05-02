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
const NUM_GENERATION = utilities.parse_int_float64(Int,readline()) 

println()

println("provide an initial temperature")
const TEMP = utilities.parse_int_float64(Float64,readline())

println()

println("provide a grid size")
const N_GRID = utilities.parse_int_float64(Int,readline())

println()

function main()

   ising_model = isingMethods.isingModel(TEMP,N_GRID) #ising model struct instantiation 
   
   isingMethods.choose_flip_strategy(ising_model)
   isingMethods.choose_trans_dynamics(ising_model)
 
   #= Global energy and magnetization time series realization will be saved on subdirectories over folder simultations=#
   mkdir("simulations/energy")
   mkdir("simulations/magnetization")

   #= Subdirectories containg a .txt file with the unicode representation of how the spin grid evolves with each generation at each run =#
   mkdir("simulations/grid_evolution")
   
   file_names = ["global_energy","global_magnetization"]

   for run in 1:NUM_RUNS      
      isingMethods.reset_stats(ising_model) 
      isingMethods.set_magnetization(INIT_MAGN,ising_model) 
      isingMethods.update_magnetization(ising_model)
      isingMethods.update_energy(ising_model)
      
      #= Generic files containing energy and magnetization time series =#
      generic_energy_file_name = "$(file_names[1])_r$(run)"*".txt"
      touch("..simulations/energy/$(generic_energy_file_name)")
      
      generic_magnetization_file_name = "$(file_names[2])_r$(run)"*".txt"
      touch("..simulations/magnetization/$(generic_magnetization_file_name)")
      
      generic_spin_grid_file_name = "spin_evolution_r$(run)"*".txt"
      touch("../simulations/$(generic_spin_grid_file_name)")

      #= Initial observations of the global energy and and magnetizaton are saved to their respective .txt files=#
      generic_energy_file = open(generic_energy_file_name,"w+")
      write(generic_energy_file,"gen $(ising_model.cur_gen) \n")
      write(generic_energy_file_name, "$(ising_model.global_energy)")
      close(generic_energy_file)

      generic_magnetization_file = open(generic_magnetization_file_name,"w+")
      write(generic_energy_file,"gen $(ising_model.cur_gen) \n")
      write(generic_energy_file_name, "$(ising_model.global_energy)")
      close(generic_energy_file)
      

      for generation in 1:NUM_GENS
         isingMethods.do_generation(ising_model)
         setfield!(ising_model,:cur_gen, generation)

         #= Writing observation i to the global energy  and magnetization time series=# 
         generic_energy_file = open(generic_energy_file_name,"w+")
         write(generic_energy_file,"gen $(ising_model.cur_gen) \n")
         write(generic_energy_file_name, "$(ising_model.global_energy)")
         close(generic_energy_file)

         generic_magnetization_file = open(generic_magnetization_file_name,"w+")
         write(generic_magnetization_file,"gen $(ising_model.cur_gen) \n")
         write(generic_magnetization_file,"$(ising_model.global_magnetization)")
         close(generic_magnetization_file)
      end

      
   end  
end

main()





