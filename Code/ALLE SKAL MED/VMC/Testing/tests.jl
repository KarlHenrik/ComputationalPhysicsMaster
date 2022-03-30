include("../initializer.jl")
include("../Wavefunctions/wavefunctions.jl")
include("../hamiltonians.jl")
include("../sampler.jl")
include("../metropolis.jl")

using Test
import Random
import StaticArrays as sa
import LinearAlgebra as la

@testset "SimpleGaussian" begin
    positions = [0.0  1.0]
    wf = SimpleGaussian(0.5, [1, 1, 1])
    @test kinetic(positions, wf, temp_arr) == 
end


@testset "initialize" begin
    rng = Random.MersenneTwister()
    
    positions = initialize(10, 3, rng)
    @test size(positions) == (10, 3)
    
    positions = initialize(2, 7, rng)
    @test size(positions) == (2, 7)
    
    positions = initialize(1, 11, rng)
    @test size(positions) == (1, 11)
    
    dim = 3
    num = 15
    wf = Correlated(0.5, 0.433, [1.0, 1.0, 1.0])
    positions = initialize(dim, num, rng, wf)
    particlesNotTooClose = true
    for p1 in 1:num
        for p2 in 1:num
            if p1 != p2
                if la.norm(positions[:, p1] - positions[:, p2]) < wf.a
                    particlesNotTooClose = false
                end
            end
        end
    end
    @test particlesNotTooClose == true
end