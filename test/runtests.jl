using TopologicalNumbers
using LinearAlgebra
using Test

using CairoMakie
# using GLMakie
using Aqua

Aqua.test_all(TopologicalNumbers; ambiguities=false)

@testset "TopologicalNumbers.jl" begin
    @testset "1D case" begin
        function H₀(k, p)

            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
        H(k) = H₀(k, (0.9, 1.0))

        N = 51
        k = range(-π, π, length=N)
        bandsum = (-60.89985395515161, 60.89985395515161)
        result = showBand(H)

        @test result.k == k
        for i in 1:size(result.Ene, 2)
            @test sum(result.Ene[:, i]) ≈ bandsum[i]
        end

        @test abs(sum(result.Ene)) < 1e-10

        @test calcBerryPhase(H) == (TopologicalNumber=[1, 1], Total=0)

        @test norm(calcBerryPhase(H).TopologicalNumber - calcBerryPhase(H, rounds=false).TopologicalNumber) < 1e-10


        H(k, p) = H₀(k, (p, 1.0))

        param = range(-2.0, 2.0, length=11)
        result = calcPhaseDiagram(H, param, "BerryPhase")

        num = [0 0; 0 0; 0 0; 1 1; 1 1; 1 1; 1 1; 1 1; 0 0; 0 0; 0 0]
        @test result.nums == num

        fig = plot1D(result.nums, result.param; disp=false)
        @test typeof(fig) == Makie.Figure

        fig = plot1D(result; disp=false)
        @test typeof(fig) == Makie.Figure

        result = calcPhaseDiagram(H, param, "BerryPhase"; rounds=false)
        num = [2.8271597168564594e-16 1.987846675914698e-17; 0.0 1.766974823035287e-17; 5.654319433712919e-16 3.533949646070574e-17; 0.9999999999999989 1.0000000000000002; 1.000000000000001 1.0000000000000002; 1.0000000000000018 1.0; 1.0 1.0000000000000007; 1.0 0.9999999999999999; 5.654319433712919e-16 7.50964299789997e-17; 2.8271597168564594e-16 1.987846675914698e-17; 5.654319433712919e-16 3.3130777931911632e-18]
        @test result.nums ≈ num


        param = range(-2.0, 2.0, length=3)
        result = calcPhaseDiagram(H₀, param, param, "BerryPhase")
        num = zeros(2, 3, 3)
        num[:, :, 1] = [0 1 0; 0 1 1]
        num[:, :, 2] = [0 0 0; 0 0 0]
        num[:, :, 3] = [0 1 0; 1 1 0]
        # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
        @test result.nums == num

        fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
        @test typeof(fig) == Makie.Figure

        fig = plot2D(result; disp=false)
        @test typeof(fig) == Makie.Figure

        param = range(-2.0, 2.0, length=4)
        result = calcPhaseDiagram(H₀, param, param, "BerryPhase"; rounds=false)
        num = zeros(2, 4, 4) # ↓結果おかしい気がする
        num[:, :, 1] = [8.481479150569378e-16 0.9999999999999994 0.9999999999999994 0.4803921568627457; 2.760898160992636e-16 0.9999999999999999 0.9999999999999997 0.9999999999999994]
        num[:, :, 2] = [0.0 5.654319433712919e-16 0.4803921568627457 0.0; 2.2639364920139617e-17 2.7056801977727833e-16 0.9999999999999994 1.6565388965955818e-17]
        num[:, :, 3] = [0.0 0.4803921568627457 5.654319433712919e-16 0.0; 1.6565388965955818e-17 0.9999999999999994 2.7056801977727833e-16 2.2639364920139617e-17]
        num[:, :, 4] = [0.4803921568627457 0.9999999999999994 0.9999999999999994 8.481479150569378e-16; 0.9999999999999994 0.9999999999999997 0.9999999999999999 2.760898160992636e-16]
        @test result.nums ≈ num

    end

    @testset "2D case" begin

        @testset "Chern" begin
            @testset "Square Flux" begin
                function H₀(k, p) # landau
                    k1, k2 = k
                    t = 1

                    Hsize = 6
                    Hmat = zeros(ComplexF64, Hsize, Hsize)

                    ϕ = 2π * p / Hsize

                    for i in 1:Hsize
                        Hmat[i, i] = -2t * cos(k2 - i * ϕ)
                    end

                    for i in 1:Hsize-1
                        Hmat[i, i+1] = -t
                        Hmat[i+1, i] = -t
                    end

                    Hmat[1, Hsize] = -t * exp(-im * Hsize * k1)
                    Hmat[Hsize, 1] = -t * exp(im * Hsize * k1)

                    Hmat
                end
                H(k) = H₀(k, 1)

                N = 51
                k = range(-π, π, length=N)
                bandsum = (-8027.411637742611, -3926.813524046014, -1090.4729277337524, 1090.4729277337524, 3926.813524046014, 8027.411637742611)
                result = showBand(H)

                @test result.k[:, 1] == k
                @test result.k[:, 2] == k
                for i in 1:size(result.Ene, 3)
                    @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
                end

                @test abs(sum(result.Ene)) < 1e-10

                @test calcChern(H) == (TopologicalNumber=[6, 6, -12, -12, 6, 6], Total=0)


                # H(k, p) = H₀(k, (p, 1.0))

                param = 1:6 # range(1, 6, length=6)
                result = calcPhaseDiagram(H₀, param, "Chern")

                #↓結果おかしい
                num = [6 6 -12 -12 6 6; 0 6 -6 -6 6 0; 0 4 -2 -1 1 0; 0 -6 6 6 -6 0; -6 -6 12 12 -6 -6; 0 2 -6 1 0 0]
                @test result.nums == num

                fig = plot1D(result.nums, result.param; disp=false)
                @test typeof(fig) == Makie.Figure

                fig = plot1D(result; disp=false)
                @test typeof(fig) == Makie.Figure

                #↓結果おかしい
                result = calcPhaseDiagram(H₀, param, "Chern"; rounds=false)
                num = [6.000000000000002 6.0 -12.0 -12.0 6.0000000000000036 5.9999999999999964; -1.0784726957737412e-16 6.0 -6.0 -6.000000000000002 5.999999999999998 4.440892098500626e-16; 2.1763752123554674e-15 4.0 -2.0000000000000018 -1.000000000000003 1.000000000000001 -1.2492169914176273e-16; 5.358542476396547e-16 -6.0 5.999999999999998 6.0 -6.000000000000002 -3.3306690738754696e-15; -6.0000000000000036 -6.0 12.0 11.999999999999996 -6.0 -6.0; 4.018335873674939e-19 2.0 -5.999999999999998 0.9999999999999992 1.1113212383547472e-15 2.5759430987039272e-15]
                @test result.nums ≈ num

                # param = range(-2.0, 2.0, length=3)
                # result = calcPhaseDiagram(H₀, param, param, "Chern")
                # num = zeros(2, 3, 3)
                # num[:, :, 1] = [1 0 1; 1 0 1]
                # num[:, :, 2] = [1 0 1; 1 0 1]
                # num[:, :, 3] = [1 0 1; 1 0 1]
                # # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
                # @test result == num

                # fig = plot2D(result[1, :, :], param, param; disp=false)
                # @test typeof(fig) == Makie.Figure

                # result = calcPhaseDiagram(H₀, param, param, "Chern"; rounds=false)
                # num = zeros(2, 3, 3)
                # num[:, :, 1] = [0.9999999999999997 0.0 1.0; 0.9999999999999982 0.0 1.0000000000000009]
                # num[:, :, 2] = [0.9999999999999997 0.0 0.9999999999999997; 0.9999999999999991 0.0 0.9999999999999991]
                # num[:, :, 3] = [1.0 0.0 1.0; 1.0000000000000018 0.0 1.0000000000000018]
                # @test result ≈ num
            end

            @testset "Haldane model" begin

                function H₀(k, p) # landau
                    k1, k2 = k
                    J = 1.0
                    K = 1.0
                    ϕ, M = p

                    h0 = 2K * cos(ϕ) * (cos(k1) + cos(k2) + cos(k1 + k2))
                    hx = -J * (1 + cos(k1) + cos(k2))
                    hy = -J * (-sin(k1) + sin(k2))
                    hz = M + 2K * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))

                    s0 = [1 0; 0 1]
                    sx = [0 1; 1 0]
                    sy = [0 -im; im 0]
                    sz = [1 0; 0 -1]

                    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
                end
                H(k) = H₀(k, (π / 2, 1.0))

                N = 51
                k = range(-π, π, length=N)
                bandsum = (-7699.060670494975, 7699.060670494976)
                result = showBand(H)

                @test result.k[:, 1] == k
                @test result.k[:, 2] == k
                for i in 1:size(result.Ene, 3)
                    @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
                end

                @test abs(sum(result.Ene)) < 1e-10

                @test calcChern(H) == (TopologicalNumber=[-1, 1], Total=0)


                H(k, p) = H₀(k, (p, 2.5))

                param = range(-π, π, length=10)
                result = calcPhaseDiagram(H, param, "Chern")

                num = [0 0; 1 -1; 1 -1; 1 -1; 0 0; 0 0; -1 1; -1 1; -1 1; 0 0]
                @test result.nums == num

                fig = plot1D(result.nums, result.param; disp=false)
                @test typeof(fig) == Makie.Figure

                fig = plot1D(result; disp=false)
                @test typeof(fig) == Makie.Figure


                result = calcPhaseDiagram(H, param, "Chern"; rounds=false)
                num = [1.1109801146677412e-15 -2.850714348572913e-16; 0.9999999999999996 -0.9999999999999991; 1.0000000000000007 -1.0000000000000002; 1.0 -1.000000000000001; 2.196996748971226e-15 4.455768167494754e-18; 6.938689780587844e-16 1.9512837637595136e-16; -0.999999999999998 0.9999999999999998; -0.9999999999999998 1.0; -1.0 0.9999999999999997; 6.668909048176786e-16 -2.7518073289070415e-16]
                @test result.nums ≈ num

                param1 = range(-π, π, length=6)
                param2 = range(-6.0, 6.0, length=6)
                result = calcPhaseDiagram(H₀, param1, param2, "Chern")
                num = zeros(2, 6, 6)
                num[:, :, 1] = [0 0 0 0 0 0; 0 0 0 0 0 0]
                num[:, :, 2] = [0 1 0 0 -1 0; 0 -1 0 0 1 0]
                num[:, :, 3] = [0 1 1 -1 -1 0; 0 -1 -1 1 1 0]
                num[:, :, 4] = [0 1 1 -1 -1 0; 0 -1 -1 1 1 0]
                num[:, :, 5] = [0 1 0 0 -1 0; 0 -1 0 0 1 0]
                num[:, :, 6] = [0 0 0 0 0 0; 0 0 0 0 0 0]
                # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
                @test result.nums == num

                fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
                @test typeof(fig) == Makie.Figure

                fig = plot2D(result; disp=false)
                @test typeof(fig) == Makie.Figure

                result = calcPhaseDiagram(H₀, param1, param2, "Chern"; rounds=false)
                num = zeros(2, 6, 6)
                num[:, :, 1] = [-1.3541921117431102e-16 -2.8652787045839576e-16 -2.940035673152271e-17 -2.3234916503891985e-16 -5.732914855719844e-17 -3.685200808344823e-16; 2.3850592097202775e-16 2.9976771205590897e-15 -4.449225133094178e-17 1.7857757969342297e-16 7.122836834784753e-16 2.3850592097202775e-16]
                num[:, :, 2] = [-6.396964825861168e-16 0.9999999999999999 -4.386804391863469e-16 1.485636657747227e-16 -0.9999999999999996 -6.388498790484686e-16; 2.458841426448316e-15 -1.0 4.2851116551941513e-16 1.2963818884457802e-17 1.0000000000000018 2.0147522165982533e-15]
                num[:, :, 3] = [4.8478557962159796e-17 1.0 0.9999999999999998 -0.9999999999999982 -0.999999999999999 2.9389535520921715e-17; 6.750172203032566e-16 -1.0 -0.9999999999999998 1.000000000000002 1.0000000000000007 4.529726153782252e-16]
                num[:, :, 4] = [4.528604538904347e-16 1.0000000000000013 1.0000000000000022 -0.9999999999999989 -0.9999999999999991 6.74905058815466e-16; 2.3723648804926086e-17 -0.999999999999999 -0.9999999999999989 0.9999999999999998 1.0 5.532730494153866e-17]
                num[:, :, 5] = [2.015700196865251e-15 1.0000000000000018 1.9702135958630458e-17 1.5461001801600323e-15 -1.0000000000000009 2.4597894067153135e-15; -8.92384864767788e-16 -1.0000000000000002 1.5185712878849862e-16 -3.9554093358765464e-16 0.9999999999999997 -8.745404034506849e-16]
                num[:, :, 6] = [4.607397984077055e-16 1.5706863748485992e-16 7.33578008994055e-16 7.325091046230339e-16 2.22052100332148e-15 4.607397984077055e-16; -3.7796154778268895e-16 -5.0747412694646605e-17 -2.4612099824722765e-16 -2.1046482842831816e-17 -3.218757251068939e-16 -3.757530988728557e-16]
                @test result.nums ≈ num
            end

        end

        @testset "Z2" begin

            function H₀(k, p) # 2d Kane-Mele
                k1, k2 = k
                t, λₛ = p

                R3 = 2λₛ * (sin(k1) - sin(k2) - sin(k1 - k2))
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
            H(k) = H₀(k, (1.0, 1.0))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-7404.378662190171, -7404.378662190169, 7404.378662190169, 7404.378662190172)
            result = showBand(H)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end

            @test abs(sum(result.Ene)) < 1e-10

            @test calcZ2(H) == (TopologicalNumber=[1, 1], Total=0)

            @test norm(calcZ2(H, rounds=false).TopologicalNumber - calcZ2(H).TopologicalNumber) < 1e-10

            @test calcZ2(H, TR=true) == (TopologicalNumber=[1, 1], TRTopologicalNumber=[1, 1], Total=0)

            @test norm(calcZ2(H, rounds=false, TR=true).TRTopologicalNumber - calcZ2(H, TR=true).TRTopologicalNumber) < 1e-10


            H(k, p) = H₀(k, (p, 1.0))

            param = range(-2.0, 2.0, length=11)
            result = calcPhaseDiagram(H, param, "Z2")

            num = [1 1; 1 1; 1 1; 1 1; 1 1; 0 0; 1 1; 1 1; 1 1; 1 1; 1 1]
            @test result.nums == num

            fig = plot1D(result.nums, result.param; disp=false)
            @test typeof(fig) == Makie.Figure

            fig = plot1D(result; disp=false)
            @test typeof(fig) == Makie.Figure

            result = calcPhaseDiagram(H, param, "Z2"; rounds=false)
            num = [1.0000000000000004 1.0; 1.0 1.0; 1.0000000000000002 1.0000000000000009; 1.0 0.9999999999999991; 0.9999999999999997 1.0; 0.0 0.0; 0.9999999999999997 1.0; 1.0000000000000004 0.9999999999999991; 1.0000000000000002 1.0000000000000009; 1.0 1.0; 1.0000000000000004 1.0]
            @test result.nums ≈ num


            param = range(-2.0, 2.0, length=3)
            result = calcPhaseDiagram(H₀, param, param, "Z2")
            num = zeros(2, 3, 3)
            num[:, :, 1] = [1 0 1; 1 0 1]
            num[:, :, 2] = [1 0 1; 1 0 1]
            num[:, :, 3] = [1 0 1; 1 0 1]
            # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
            @test result.nums == num

            fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
            @test typeof(fig) == Makie.Figure

            fig = plot2D(result; disp=false)
            @test typeof(fig) == Makie.Figure

            result = calcPhaseDiagram(H₀, param, param, "Z2"; rounds=false)
            num = zeros(2, 3, 3)
            num[:, :, 1] = [0.9999999999999997 0.0 1.0; 0.9999999999999982 0.0 1.0000000000000009]
            num[:, :, 2] = [0.9999999999999997 0.0 0.9999999999999997; 0.9999999999999991 0.0 0.9999999999999991]
            num[:, :, 3] = [1.0 0.0 1.0; 1.0000000000000018 0.0 1.0000000000000018]
            @test result.nums ≈ num

        end

    end

end

