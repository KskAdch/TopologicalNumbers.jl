using TopologicalNumbers
using Test

using CairoMakie
using Aqua

# Aqua.test_all(TopologicalNumbers; ambiguities=false)

@testset "TopologicalNumbers.jl" begin
    @testset "1D case" begin
        function H(k)
            g = 0.9
            
            [
                0 g+exp(-im*k)
                g+exp(im*k) 0
            ]
        end

        fig = showBand(H, 1)
        @test typeof(fig) == Makie.Figure

        @test calcBerryPhase(H) == (TopologicalNumber = [1, 1], Total = 0)
    end

    @testset "2D case" begin

        @testset "Chern" begin
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

            fig = showBand(H, 2)
            @test typeof(fig) == Makie.Figure

            @test calcChern(H) == (TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)
        end

        @testset "Z2" begin
            function H(k) # 2d Kane-Mele
                k1, k2 = k
                t = 1

                R3 = 2(sin(k1) - sin(k2) - sin(k1-k2))
                R4 = -t*(sin(k1) + sin(k2))
                R5 = -t*(cos(k1) + cos(k2) + 1)

                s0 = [1 0; 0 1]
                sx = [0 1; 1 0]
                sy = [0 -im; im 0]
                sz = [1 0; 0 -1]

                a3 = kron(sz, sz)
                a4 = kron(sy, s0)
                a5 = kron(sx, s0)

                R3*a3+R4*a4+R5*a5
            end

            fig = showBand(H, 2)
            @test typeof(fig) == Makie.Figure

            @test calcZ2(H) == (TopologicalNumber = [1, 1], Total = 0)

            @test calcZ2(H, TR = true) == (TopologicalNumber = [1, 1], TRTopologicalNumber = [1, 1], Total = 0)
        end

    end

    # test


end
