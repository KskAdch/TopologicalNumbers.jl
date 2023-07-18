using TopologicalNumbers
using Test

@testset "TopologicalNumbers.jl" begin
    @testset "1D case" begin
        function H(k)
            g = 0.9
            
            [
                0 g+exp(-im*k)
                g+exp(im*k) 0
            ]
        end
        @test QuantizedBerryPhase(H) == (TopologicalNumber = [1, 1], Total = 0)
    end

    @testset "2D case" begin
        function H(k) # landau
            k1, k2 = k
            t = 1
        
            Hsize = 6
            Hmat = zeros(ComplexF64, Hsize, Hsize)
        
            for i in 1:Hsize
                Hmat[i, i] = -2*cos(k2-2pi*i/Hsize)
            end
        
            for i in 1:Hsize-1
                Hmat[i, i+1] = -t
                Hmat[i+1, i] = -t
            end
        
            Hmat[1, Hsize] = -t*exp(-im*k1)
            Hmat[Hsize, 1] = -t*exp(im*k1)
            
            Hmat
        end

        @test FirstChern(H) == (TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)

    end

    # test


end
