module isingMethods
export display,reset_stats,compute_energy_cell,update_energy,update_magenetization,randomize,set_magnetization,try_cell_flip,get_cell_coords,get_cell_id
export do_generation, choose_flip_strategy, choose_trans_dynamics

include("../src/ising.jl")
using .ising: isingModel,CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS,GLAUBER_DYNAMICS
export isingModel, CRITICAL_TEMP, RANDOM_STRATEGY, SHUFFLE_STRATEGY, SEQUENTIAL_STRATEGY, METROPOLIS_DYNAMICS,GLAUBER_DYNAMICS

include("utilities.jl")
using .utilities: swap!

using Random 
Random.seed!(1234)

#=Unicode representation of the spin grid=#
function display(ising_model :: isingModel)
    for i in 1:ising_model.NGRID
        for j in 1:ising_model.NGRID
            if ising_model.grid[i,j] === 1
                print("+")
            else
                print("-")
            end     
        end
        println()
    end    
end

#=Resets statistics=#
function reset_stats(ising_model :: isingModel)
    fields_names_to_reset = ["global_energy","global_mean","global_variance",
    "global_magnetization"]

     for (i,field_name_string) in enumerate(fields_names_to_reset)
        field_name = Symbol(field_name_string) #converts string to Symbol 
        setfield!(ising_model,field_name,0.0)
     end      
end


#=Returns the energy of a cell
The grid wraps around at the edges (toroidal symmetry).=#
function compute_energy_cell(i :: Int, j :: Int,ising_model :: ising.isingModel) :: Float64
    grid_spin = ising_model.grid  
    energy = 0
    ip, im = 1,1 
    jp, jm = 1,1

    #Toroidal symetry 
    if i === 1 
        im = ising_model.NGRID
    elseif i === ising_model.NGRID
        ip = 1  
    end

    if j === 1 
        jm = ising_model.NGRID
    elseif j === ising_model.NGRID
        jp = 1
    end
    
    #Point (i,j) is not on the lattice boundries 
    if (1 < i < ising_model.NGRID) && (1 < j < ising_model.NGRID)
        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
    end
    
    #Energy of cell i,j is determined by the energy of its neighbours 
    energy += grid_spin[ip,j] + grid_spin[im,j] + grid_spin[i,jp] + grid_spin[i,jm]
    energy = -grid_spin[i,j]*energy
    return energy #energy per cell 
end

#=Computes the enery of the whole grid of spins states=#
function update_energy(ising_model :: isingModel)
    g_energy = 0
    for i in 1:ising_model.NGRID
        for j in 1:ising_model.NGRID
            ising_model.global_energy += compute_energy_cell(i,j,ising_model)
            g_energy = ising_model.global_energy
        end
    end
   setfield!(ising_model, :global_energy,g_energy)
end

#=Computes the magnetization of the whole spin grid=#
function update_magenetization(ising_model :: isingModel)
    g_magnetization = 0
    for i in 1:ising_model.NGRID
        for j in 1:ising_model.NGRID
            ising_model.global_magnetization += ising_model.grid[i,j]
            g_magnetization = ising_model.global_magnetization
        end
    end

    g_magnetization /= ising_model.NGRID
    setfield!(ising_model,:global_magnetization,g_magnetization)
end

#=Randomly populates the spin grid=#
function randomize(ising_model :: isingModel)
    spin_grid = ising_model.grid
    rn = 0 
    for i in 1:ising_model.NGRID
        for j in 1:ising_model.NGRID
          rn = rand(Int)
          if rn % 2 === 0 
            spin_grid[i,j] =-1;
          else 
            spin_grid[i,j] = 1;
          end
        end
    end
    setfield!(ising_model,:grid,spin_grid)
end

#=Generates random changes in the spin grid but conserving a given magnetization=#
function set_magnetization(magn :: Float64, ising_model :: isingModel)
    spin_grid = ising_model.grid
    p = (1+magn)/2.0
    
    for i in 1:ising_model.NGRID
        for j in 1:ising_model.NGRID
            if rand() <= p 
                spin_grid[i,j] = +1
            else
                spin_grid[i,j] = -1   
            end
        end
    end
    setfield!(ising_model,:grid,spin_grid)
end

#=Method for trying a flip spin=#
function try_cell_flip(i :: Int, j :: Int, ising_model ::isingModel)
    g_energy = 0
    prob = 0
    old_E = compute_energy_cell(i,j,ising_model)
    new_E = -old_E #Always true since E_i = s_i*(sum_neighs s_n)
    ΔE = new_E - old_E 
    if ising_model.trans_dynamics === ising.METROPOLIS_DYNAMICS
        #Metropolis dynamics
        if ΔE <= 0
            prob = 1 #energy is lower, spin is flipped 
        else
            prob = exp(-ΔE/ising_model.TEMP) #termal flip 
        end
    elseif ising_model.trans_dynamics === ising.DYNAMICS_GLAUBER
        #glauber dynamics
        prob = 1/(1 + exp(-ΔE/ising_model.TEMP))   
    end  
    
    if rand() <= prob
        do_flip = true 
    else
        do_flip = false    
    end
        
    if do_flip
        #determine changes in global magnetization
        if ising_model.grid[i,j] === 1
            ising_model.global_magnetization -= 2/ising_model.NCELLS
        else
            ising_model.global_magnetization += 2/ising_model.NCELLS    
        end 

        ising_model.grid[i,j] *=-1 #cell gets flipped
        spin_grid = ising_model.grid
        setfield!(ising_model,:grid,spin_grid) 

        ising_model.global_energy += ΔE #individual energy changed
        g_energy = ising_model.global_energy
        setfield!(ising_model,:global_energy,g_energy) 
    end     

end

#=
NOTES: 
1) Because Julia is 1-indexed in the folowing the base function mod1(x,y) will be used. This function returns an integer r  
in the set (0,y] i.e is the same as the % operator but with an offset of 1

2) the ceiled divison function cld(x,y) will be used. 
=#

#= Given the id that uniquely determines a cell in the spin grid, the function outputs the (i,j) 
coordinates inside the grid location=#
function get_cell_coords( id :: Int,ising_model :: isingModel) :: Array{Int,1}
    i = mod1(id ,ising_model.NGRID)
    j = cld(id,ising_model.NGRID)
    return [i,j]
end

#=Provided the (x,y) coordinates of a cell gives the id representation of a spin at location (x,y)=#
function get_cell_id(i :: Int, j :: Int, ising_model :: isingModel) :: Int
    return i + (j-1)*ising_model.NGRID   
end

#=Applies a flip cell spins to each spin in teh spin grid=#
function do_generation(ising_model :: isingModel)
    # with strategy it may happen that not all cells get flipped 
    if ising_model.flip_strategy === ising.RANDOM_STRATEGY
        for temp in 1:ising_model.NCELLS #this loops from 1 to N² (the number of cells )
            i = mod1(rand(Int),ising_model.NGRID)
            j = mod1(rand(Int),ising_model.NGRID)
            try_cell_flip(i,j,ising_model) 
        end
    #all cells are granted to be flipped at least once (Fisher-Yates algorithm)
    elseif ising_model.flip_strategy === ising.SHUFFLE_STRATEGY
        ising_model.flip_order = 1:(ising_model.NGRID*ising_model.NGRID)
        fliping_order = ising_model.flip_order
        for temp in 1:ising_model.NCELLS
            i = mod1(rand(Int), ising_model.NCELLS - temp + 1 ) + temp - 1 
            utilities.swap!(temp,i,fliping_order)
        end
        
        setfield!(ising_model,:flip_order,fliping_order)

        #all cells are flipped 
        for i in 1:ising_model.NCELLS
            array_coords = get_cell_coords(ising_model.flip_order[i],ising_model)
            try_cell_flip(array_coords[1], array_coords[2], ising_model)
        end
    
    elseif ising_model.flip_strategy === ising.SEQUENTIAL_STRATEGY
        for i in 1:ising_model.NGRID
            for j in 1:ising_model.NGRID
                try_cell_flip(i,j,ising_model) #sequential flip
            end        
        end
    end
end

#=Method which allows the user to select the a flip strategy=#
function choose_flip_strategy(ising_model :: isingModel)
    options_to_choose = Dict("1" => ising.RANDOM_STRATEGY, "2" => ising.SHUFFLE_STRATEGY,
    "3" => ising.SEQUENTIAL_STRATEGY)

    while true 
        println("Choose 1 of the 3 possible flip strategies")
        println("1. Random strategy")
        println("2. Shuffle strategy")
        println("3. Sequential strategy")
        println()

        user_input = readline()
        user_choice = parse(Int64, user_input)
        user_choice_key = string(user_choice)
        
        if 1 <= user_choice <= 3
            setfield!(ising_model,:flip_strategy,get(options_to_choose,user_choice_key,nothing))
            break 
        end
    end
end

#=Method for selecting a transition dynamics=#
function choose_trans_dynamics(ising_model :: isingModel)
    options_to_choose = Dict("1" => METROPOLIS_DYNAMICS, "2" => ising.GLAUBER_DYNAMICS)

    while true 
        println("Choose a transition dynamics")
        println("1. Metropolis dynamics")
        println("2. Glauber dynamics")
        println()

        user_input = readline()
        user_choice = parse(Int64, user_input)
        user_choice_key = string(user_choice)
        
        if 1 <= user_choice <= 2
            setfield!(ising_model,:trans_dynamics,get(options_to_choose,user_choice_key,nothing))
            break 
        end
    end
end

end 

