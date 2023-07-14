#= This will test the numerical convergence of the simulations at temperatures below, far and equal to the
critical temperature =#     

using Test 

include("../scripts/src/utilities.jl")
using utilities: get_array_from_txt


magnTBelowTc = [0.387712431593414,0.359732610984483,0.35397601016284,0.341585752687131,0.39589851829128,0.364744771633911,
0.376254820064415,0.323906345916509,0.323826618097163,0.327263782748425,
0.239645612244689,0.181323571532698,0.107096621179225] #array with magnetization at temperatures in the interval 1 ≤ T ≤ 2.2

magnAtTc = 0.118454602928595 # T = Tc 

magnFarTc = [0.168990575856311,0.0589899351034247,0.0207203717088545,0.00994495400708741,0.00665357733560256,
0.00590289351320215,0.00386716123174351,0.0031434900327604,0.00291024033488576,0.00226280509445688,0.00188265867306743
,0.00180068581052776,0.00131281131627664] #Tc < T ≤ 3 

@testset "magnetization at temperatures below Tc" begin
    @test 
    @test 
    @test 
    @test 
    @test 
    @test 
    @test 
    @test 
    @test
    @test
    @test
    @test
    @test
end

@testset "magnetization at Tc" begin
    @test 
end

@testset "magnetization far Tc" begin
    @test 
    @test 
    @test 
    @test 
    @test 
    @test 
    @test 
    @test 
    @test
    @test
    @test
    @test
    @test
end
