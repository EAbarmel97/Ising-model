#Needed imports 
include("utilities.jl"
,"flippingStrategies.jl")


#Parameters needed to construct the ising model 
mutable struct isingModel
    #The critical temperature is Tc = 2/ln(1+sqrt(2)) in units of J/k, with k beeing the Boltzman constant 
    const CRITICAL_TEMP = 2.26918531421302

    #model constants 
    TEMP::Double64
    NGRID::Int
    NCELL ::Int

    grid :: Array{Int,2} #2D array which will contains the spin vals 
    gridCopy :: Array{Int,2} #copy of the spin 2D array 
    
    #Dead cells and density of dead cells 
    deadCell :: Array{Bool,2}
    useDeadCells :: Bool
    deadDensity :: Double64

    #Fliping strategies 
    flipStrategy:: Int
    const RANDOM_STRATEGY = 0
    const SHUFFLE_STRATEGY = 1
    const SEQUENTIAL_STRATEGY = 2

    #Transition dynamics 
    const METROPOLIS_DYNAMICS = 0

    #Global statistics
    globalEnergy :: Double64
    globalMean :: Double64
    gobalVariance :: Double64
    globalMagnetization :: Double64
    globalM2 :: Double64 
    globalNpoints :: Int 

    #Tracking samples 
    trackSamples :: Bool

    #Sampling variables and sampling constants 
    NUM_SAMPLES :: Int
    SAMPLE_MIN :: Int 
    SAMPLE_MAX :: Int

    sampleSize :: Int
    #optional sampling variables 
    sampleMagnetization :: Double64
    sampleMean :: Double64
    sample_var :: Double64
    sampleM2 :: Double64 
    sampleNpts :: Double64
    
    #number of points to remember 
    NUM_DATA :: Int 
    
    #Running var/ mean 
    runningMean :: Double64
    runningVar :: Double64

    nextData :: Int 

    #model constructor with two params 
    function isingModel(pTemp)

    













        return new(

        )    
    end

    function ()
        
    end
    
end

function getCellCoords!(ID:: Int, x::Int, y::Int)
    
end

function updateSampleMagnetization(sample :: Int)
    sum = 0 
    for i in 1:sampleSize[sample]
        getCellCoords()
end
