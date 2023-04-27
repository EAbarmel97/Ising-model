using CSV, Tables #libraries to export arrays and save them as CSV files

include("../scripts/.julia/src/isingMethods.jl")
using .isingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS,GLAUBER_DYNAMICS
using .isingMethods: display,reset_stats,compute_energy_cell,update_energy,update_magenetization,randomize,set_magnetization
using .isingMethods: get_cell_coords, get_cell_id, do_generation,choose_flip_strategy

include("../scripts/.julia/src/utilities.jl")
using .utilities: parse_int_float64

println("provide an initial magnetization")
const INIT_MAGN = utilities.parse_int_float64(Float64,readline())
 
println("provide the number of runs")
const NUM_RUNS = utilities.parse_int_float64(Int,readline())

println("provide the number of simulations")
const NUM_SIMULATIONS = utilities.parse_int_float64(Int,readline()) 

println("provide an initial temperature")
const TEMP = utilities.parse_int_float64(Float64,readline())

println("provide a grid size")
const N_GRID = utilities.parse_int_float64(Int,readline())

function main()

   ising_model = isingMethods.isingModel(TEMP,N_GRID) #ising model struct instantiation 
   
   isingMethods.choose_flip_strategy(ising_model)
   isingMethods.choose_trans_dynamics(ising_model)

   file_names = ["global_magnetization.csv","global_energy.csv"]

   #creation of CSV files that will contain the magnetization and global energy time series, respectively 
   touch("global_magnetization.csv")
   global_magnetization_csv = open("global_magnetization.csv","w+")

   touch("global_energy.csv")
   global_energy_csv = open("global_energy.csv","w+")
   

   #arrays containing the global energy and global magnetization time series
   global_energy_time_series = []
   global_magnetization_time_series = []

   
   for run in 1:NUM_RUNS
      isingMethods.reset_stats(ising_model)
      isingMethods.update_energy(ising_model)
      isingMethods.set_magnetization(INIT_MAGN,ising_model)
      isingMethods.update_magenetization(ising_model)

      for simulation in 1:NUM_SIMULATIONS
         isingMethods.do_generation(ising_model)

         #pushing observations of the time series at time i to time series arrays 
         push!(global_energy_time_series,ising_model.global_energy)
         push!(global_magnetization_time_series,ising_model.global_magnetization)

         for i in eachindex(file_names)
            file_name = file_names[i]
            if file_name === "global_magnetization.csv"
               println("Recording time series on file $file_name")
               CSV.write(file_name,Tables.table(
                  global_energy_time_series),writeheader=false)
            elseif file_name === "global_energy.csv"
               println("Recording time series on file $file_name") 
               CSV.write(file_name,Tables.table(
                  global_magnetization_time_series),writeheader=false)
            end      
         end
      end 

      close(global_magnetization_csv)
      close(global_energy_csv)
   end  
end

main()





