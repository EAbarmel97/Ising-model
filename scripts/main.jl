using CSV 

include("../scripts/.julia/src/isingMethods.jl")
using .isingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS,GLAUBER_DYNAMICS
using .isingMethods: display,reset_stats,compute_energy_cell,update_energy,update_magenetization,randomize,set_magnetization
using .isingMethods: get_cell_coords, get_cell_id, do_generation,choose_flip_strategy

include("../scripts/.julia/src/utilities.jl")
using .utilities: parse_int_float64

function main()

   println("provide a initial magnetization")
   INIT_MAGN = utilities.parse_int_float64(Float64,readline())
    
   println("provide the number of runs")
   NUM_RUNS = utilities.parse_int_float64(Int,readline())
   
   println("provide the number of simulations")
   NUM_SIMULATIONS = utilities.parse_int_float64(Int,readline()) 

   println("provide an initial temperature")
   TEMP = utilities.parse_int_float64(Float64,readline())

   println("provide a grid size")
   SIZE_GRID = utilities.parse_int_float64(Int,readline())

   ising_model = isingMethods.isingModel(TEMP,SIZE_GRID) #ising model struct instantiation 
   
   isingMethods.choose_flip_strategy(ising_model)
   isingMethods.choose_trans_dynamics(ising_model)

   #creation of CSV files that will contain the magnetization and global energy time series, respectively 
   touch("magn.csv")
   touch("global_energy.csv")

   file_names = ["magn.csv","global_energy"]
   file_name = ""
   
   for run in 1:NUM_RUNS
      isingMethods.reset_stats(ising_model)
      isingMethods.update_energy(ising_model)
      isingMethods.set_magnetization(INIT_MAGN,ising_model)
      isingMethods.update_magenetization(ising_model)

      for simulation in 1:NUM_SIMULATIONS
         isingMethods.do_generation(ising_model)

         for i in eachindex(file_names)
            file_name = file_names[i]
            println("Recording on file $file_name")
            open(file_name,"wr")

         end
         
      end   
   end
  
end

main()





