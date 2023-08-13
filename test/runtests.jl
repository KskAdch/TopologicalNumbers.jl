using TopologicalNumbers
using LinearAlgebra
using Test

# using CairoMakie
# using GLMakie
using Aqua

Aqua.test_all(TopologicalNumbers; ambiguities=false)

@testset "TopologicalNumbers.jl" begin
    @testset "1D case" begin
        function H(k)
            g = 0.9

            [
                0 g+exp(-im * k)
                g+exp(im * k) 0
            ]
        end

        N = 51
        # k = range(-π, π, length=N)
        bandsum = (-60.89985395515161, 60.89985395515161)
        result = showBand(H)

        # @test result.k .== k
        for i in 1:size(result.Ene, 2)
            @test sum(result.Ene[:, i]) ≈ bandsum[i]
        end

        @test abs(sum(result.Ene)) < 1e-10

        @test calcBerryPhase(H) == (TopologicalNumber=[1, 1], Total=0)

        @test norm(calcBerryPhase(H).TopologicalNumber - calcBerryPhase(H, rounds=false).TopologicalNumber) < 1e-10


        function H0(k, p)

            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
        H(k, p) = H0(k, (p, 1.0))

        param = range(-2.0, 2.0, length=11)
        result = calcPhaseDiagram(H, param, "BerryPhase")

        num = [0 0; 0 0; 0 0; 1 1; 1 1; 1 1; 1 1; 1 1; 0 0; 0 0; 0 0]
        @test result == num

        param = range(-2.0, 2.0, length=3)
        result = calcPhaseDiagram(H0, param, param, "BerryPhase")
        num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
        @test result == num

    end

    @testset "2D case" begin

        @testset "Chern" begin
            function H(k) # landau
                k1, k2 = k
                t = 1

                Hsize = 6
                Hmat = zeros(ComplexF64, Hsize, Hsize)

                for i in 1:Hsize
                    Hmat[i, i] = -2 * cos(k2 - 2pi * i / Hsize)
                end

                for i in 1:Hsize-1
                    Hmat[i, i+1] = -t
                    Hmat[i+1, i] = -t
                end

                Hmat[1, Hsize] = -t * exp(-im * k1)
                Hmat[Hsize, 1] = -t * exp(im * k1)

                Hmat
            end

            N = 51
            # k = range(-π, π, length=N)
            bandsum = (-8026.922381020278, -3931.3204155463714, -1076.7367579094316, 1076.7367579094316, 3931.320415546372, 8026.922381020279)
            result = showBand(H)

            # @test result.k[:, 1] .== k
            # @test result.k[:, 2] .== k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end

            @test abs(sum(result.Ene)) < 1e-10

            @test calcChern(H) == (TopologicalNumber=[1, 1, -2, -2, 1, 1], Total=0)


            # function H(k, p)
            #     k1, k2 = k
            #     t = p

            #     Hsize = 6
            #     Hmat = zeros(ComplexF64, Hsize, Hsize)

            #     for i in 1:Hsize
            #         Hmat[i, i] = -2 * cos(k2 - 2pi * i / Hsize)
            #     end

            #     for i in 1:Hsize-1
            #         Hmat[i, i+1] = -t
            #         Hmat[i+1, i] = -t
            #     end

            #     Hmat[1, Hsize] = -t * exp(-im * k1)
            #     Hmat[Hsize, 1] = -t * exp(im * k1)

            #     Hmat
            # end
            # # H(k, p) = H0(k, (p, 1.0))

            # param = range(-2.0, 2.0, length=11)
            # result = calcPhaseDiagram(H0, param, "Chern")

            # num = [0 0 0 0 1 1 1 0 0 0 0; 0 0 0 0 1 1 1 0 0 0 0]
            # @test result == num

            # param = range(-2.0, 2.0, length=3)
            # result = calcPhaseDiagram(H, param, param, "Chern")
            # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
            # @test result == num


        end

        @testset "Z2" begin
            function H(k) # 2d Kane-Mele
                k1, k2 = k
                t = 1

                R3 = 2(sin(k1) - sin(k2) - sin(k1 - k2))
                R4 = -t * (sin(k1) + sin(k2))
                R5 = -t * (cos(k1) + cos(k2) + 1)

                s0 = [1 0; 0 1]
                sx = [0 1; 1 0]
                sy = [0 -im; im 0]
                sz = [1 0; 0 -1]

                a3 = kron(sz, sz)
                a4 = kron(sy, s0)
                a5 = kron(sx, s0)

                R3 * a3 + R4 * a4 + R5 * a5
            end

            N = 51
            # k = range(-π, π, length=N)
            bandsum = (-7404.378662190171, -7404.378662190169, 7404.378662190169, 7404.378662190172)
            result = showBand(H)

            # @test result.k[:, 1] .== k
            # @test result.k[:, 2] .== k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end

            @test abs(sum(result.Ene)) < 1e-10

            @test calcZ2(H) == (TopologicalNumber=[1, 1], Total=0)

            @test norm(calcZ2(H, rounds=false).TopologicalNumber - calcZ2(H).TopologicalNumber) < 1e-10

            @test calcZ2(H, TR=true) == (TopologicalNumber=[1, 1], TRTopologicalNumber=[1, 1], Total=0)

            @test norm(calcZ2(H, rounds=false, TR=true).TRTopologicalNumber - calcZ2(H, TR=true).TRTopologicalNumber) < 1e-10


            # function H(k, p) # 2d Kane-Mele
            #     k1, k2 = k
            #     t = p

            #     R3 = 2(sin(k1) - sin(k2) - sin(k1 - k2))
            #     R4 = -t * (sin(k1) + sin(k2))
            #     R5 = -t * (cos(k1) + cos(k2) + 1)

            #     s0 = [1 0; 0 1]
            #     sx = [0 1; 1 0]
            #     sy = [0 -im; im 0]
            #     sz = [1 0; 0 -1]

            #     a3 = kron(sz, sz)
            #     a4 = kron(sy, s0)
            #     a5 = kron(sx, s0)

            #     R3 * a3 + R4 * a4 + R5 * a5
            # end
            # # H(k, p) = H0(k, (p, 1.0))

            # param = range(-2.0, 2.0, length=11)
            # result = calcPhaseDiagram(H, param, "Chern")

            # num = [0 0 0 0 1 1 1 0 0 0 0; 0 0 0 0 1 1 1 0 0 0 0]
            # @test result == num

            # param = range(-2.0, 2.0, length=3)
            # result = calcPhaseDiagram(H0, param, param, "Chern")
            # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
            # @test result == num

        end

    end

end
