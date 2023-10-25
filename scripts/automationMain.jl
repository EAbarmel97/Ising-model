using Random
Random.seed!(1234)

include("../scripts/src/IsingMethods.jl")
using .IsingMethods: isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS, GLAUBER_DYNAMICS
using .IsingMethods: display, reset_stats, compute_energy_cell, update_energy, update_magnetization, randomize, set_magnetization
using .IsingMethods: get_cell_coords, get_cell_id, do_generation, choose_flip_strategy

include("src/utils/utilities.jl")
using .utilities: parse_int_float64, get_array_from_txt

include("../scripts/src/Exceptions.jl")
using .Exceptions: IlegalChoiceException


const SIMULS_DIR = "all_simulations"
const AUTOMATED_SIMULS_DIR = joinpath(SIMULS_DIR, "automated")

if !isdir(AUTOMATED_SIMULS_DIR)
   mkpath(AUTOMATED_SIMULS_DIR)
end

#Initial temperature
const INIT_TEMP = utilities.parse_int_float64(Float64, ARGS[1])

#Final temperature
const FINAL_TEMP = utilities.parse_int_float64(Float64, ARGS[2])

if FINAL_TEMP < INIT_TEMP 
   throw(Exceptions.IlegalChoiceException("Ilegal  choice. Tf < Ti"))   
end

#Number of runs
const NUM_RUNS = utilities.parse_int_float64(Int, ARGS[3])

#Increment
const INCREMENT = utilities.parse_int_float64(Float64, ARGS[4])

#= number of different temperatures equaly spaced by given incrments, contained in the interval [Ti, Tf] =#
const NUM_TEMPS = ceil(Int,(FINAL_TEMP - INIT_TEMP)/INCREMENT)

#Grid size
const N_GRID = utilities.parse_int_float64(Int, ARGS[5])

#Number of generations
const NUM_GENERATIONS = utilities.parse_int_float64(Int, ARGS[6])

cd("all_simulations/automated") #going up in the working directory

const CURR_DIR = pwd()
println()
println("saving simulations under dir: $CURR_DIR")

function do_model(INIT_MAGN, TEMP, N_GRID)
   ising_model = isingMethods.isingModel(TEMP, N_GRID) #ising model struct instantiation
   ising_model.flip_strategy = isingMethods.RANDOM_STRATEGY
   ising_model.trans_dynamics = isingMethods.METROPOLIS_DYNAMICS
   
   ROUNDED_TEMP = round(TEMP, digits=2)
   str_temp = replace("$(ROUNDED_TEMP)", "." => "_") #stringified temperature with "." replaced by "_"
   aux_dir =  CURR_DIR * "/simulations_T_" * str_temp #folder containing simulations al temp str_temp 
   FOURIER_AUTOMATED_DIR = joinpath(aux_dir, "fourier")
   mkpath(aux_dir) #creates simulation folder
   mkpath(FOURIER_AUTOMATED_DIR)  

   #= Global magnetization time series realization will be saved on subdirectories over folder simultations=#
   global_magnetization_aux_dir = aux_dir * "/magnetization"
   mkpath(global_magnetization_aux_dir)

   for run in 1:NUM_RUNS
      isingMethods.reset_stats(ising_model)
      isingMethods.set_magnetization(INIT_MAGN, ising_model) #populates the spin grid with a given initial magnetization 
      isingMethods.update_magnetization(ising_model) #updates global magnetization 
      isingMethods.update_energy(ising_model) #updates global energy

      #= Creation of generic .txt files containing global magnetization time series =#
      generic_magnetization_file_name = global_magnetization_aux_dir * "/global_magnetization_r$(run)" * ".txt"
      touch("$generic_magnetization_file_name")

      #= Initial observations of the global magnetizaton are saved to their respective .txt files=#
      generic_magnetization_file = open(generic_magnetization_file_name, "w+")
      write(generic_magnetization_file, "$(ising_model.global_magnetization)\n") #initial observation
      close(generic_magnetization_file)

      for generation in 1:NUM_GENERATIONS
         isingMethods.do_generation(ising_model)
         setfield!(ising_model, :cur_gen, generation)

         generic_magnetization_file = open(generic_magnetization_file_name, "a+")
         write(generic_magnetization_file, "$(ising_model.global_magnetization)\n") #global magnetization observation at generation i 
         close(generic_magnetization_file)

         if generation == NUM_GENERATIONS
            isingMethods.do_generation(ising_model)
            setfield!(ising_model, :cur_gen, generation)

            generic_magnetization_file = open(generic_magnetization_file_name, "a+")
            write(generic_magnetization_file, "$(ising_model.global_magnetization)") #global magnetization observation at generation i 
            close(generic_magnetization_file)
         end
      end      
   end
end

#= do_model function wrapper to make simulations of the ising model at different temperatures =#
function do_simulations(num_temps :: Int, init_temp, final_temp)
   for i in 0:num_temps
      #= random initial temperature on the interval [-1 ,1] =#
      rand_magn = rand()*2 - 1 
      
      #= temperature increments in arithmetic porgression  =#
      temp = init_temp + (i/num_temps)*(final_temp - init_temp)

      do_model(rand_magn, temp, N_GRID) 
   end   
end

function main()
   do_simulations(NUM_TEMPS, INIT_TEMP, FINAL_TEMP) 
end

main()