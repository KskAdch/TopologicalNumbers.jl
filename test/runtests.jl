using TopologicalNumbers
using LinearAlgebra
using Test

using MPI

using PythonPlot

# using Aqua
# Aqua.test_all(TopologicalNumbers; ambiguities=false)

@testset "TopologicalNumbers.jl" begin
    @testset "1D case" begin
        function H₀(k, p)

            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
        @test H₀(0.0, (1.0, 1.0)) == SSH(0.0, (1.0, 1.0))
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

        prob = BPProblem(H)
        @test solve(prob).TopologicalNumber == [1, 1]

        H(k, p) = H₀(k, (p, 1.0))

        param = range(-2.0, 2.0, length=11)
        result = calcPhaseDiagram(H, param, "BerryPhase")
        calcPhaseDiagram(H, param, "BerryPhase"; progress=true)
        result_MPI = calcPhaseDiagram(H, param, "BerryPhase"; parallel=UseMPI(MPI))
        calcPhaseDiagram(H, param, "BerryPhase"; parallel=UseMPI(MPI), progress=true)
        @test result.nums == result_MPI.nums

        num = [0 0; 0 0; 0 0; 1 1; 1 1; 1 1; 1 1; 1 1; 0 0; 0 0; 0 0]
        @test result.nums == num

        prob = BPProblem(H)
        @test calcPhaseDiagram(prob, param).nums == num

        fig = plot1D(result.nums, result.param; disp=false)
        @test typeof(fig) == Figure
        plotclose()

        fig = plot1D(result; disp=false)
        @test typeof(fig) == Figure
        plotclose()

        result = calcPhaseDiagram(H, param, "BerryPhase"; rounds=false)
        num = [2.8271597168564594e-16 1.987846675914698e-17; 0.0 1.766974823035287e-17; 5.654319433712919e-16 3.533949646070574e-17; 0.9999999999999989 1.0000000000000002; 1.000000000000001 1.0000000000000002; 1.0000000000000018 1.0; 1.0 1.0000000000000007; 1.0 0.9999999999999999; 5.654319433712919e-16 7.50964299789997e-17; 2.8271597168564594e-16 1.987846675914698e-17; 5.654319433712919e-16 3.3130777931911632e-18]
        @test result.nums ≈ num


        param = range(-2.0, 2.0, length=3)
        result = calcPhaseDiagram(H₀, param, param, "BerryPhase")
        calcPhaseDiagram(H₀, param, param, "BerryPhase"; progress=true)
        result_MPI = calcPhaseDiagram(H₀, param, param, "BerryPhase"; parallel=UseMPI(MPI))
        calcPhaseDiagram(H₀, param, param, "BerryPhase"; parallel=UseMPI(MPI), progress=true)
        @test result.nums == result_MPI.nums

        num = zeros(2, 3, 3)
        num[:, :, 1] = [0 1 0; 0 1 0]
        num[:, :, 2] = [0 0 0; 0 0 0]
        num[:, :, 3] = [0 1 0; 0 1 0]
        # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
        @test result.nums == num

        prob = BPProblem(H₀)
        @test calcPhaseDiagram(prob, param, param).nums == num

        fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
        @test typeof(fig) == Figure
        plotclose()

        fig = plot2D(result; disp=false)
        @test typeof(fig) == Figure
        plotclose()

        param = range(-2.0, 2.0, length=4)
        result = calcPhaseDiagram(H₀, param, param, "BerryPhase"; rounds=false)
        num = zeros(2, 4, 4)
        num[:, :, 1] = [3.3306690738754696e-16 0.9999999999999998 0.9999999999999997 6.661338147750939e-16; 3.3306690738754696e-16 0.9999999999999997 0.9999999999999998 6.661338147750939e-16]
        num[:, :, 2] = [0.0 0.0 8.881784197001252e-16 0.0; 0.0 0.0 8.881784197001252e-16 0.0]
        num[:, :, 3] = [0.0 8.881784197001252e-16 0.0 0.0; 0.0 8.881784197001252e-16 0.0 0.0]
        num[:, :, 4] = [6.661338147750939e-16 0.9999999999999997 0.9999999999999998 3.3306690738754696e-16; 6.661338147750939e-16 0.9999999999999998 0.9999999999999997 3.3306690738754696e-16]
        @test result.nums ≈ num

    end

    @testset "2D case" begin
        @testset "Chern" begin
            @testset "Square Flux" begin
                function H₀(k, p) # landau
                    k1, k2 = k
                    t = 1
                    Hsize, ν = p

                    # Hsize = 6
                    Hmat = zeros(ComplexF64, Hsize, Hsize)

                    ϕ = 2π * ν / Hsize

                    for i in 1:Hsize
                        Hmat[i, i] = -2t * cos(k2 - i * ϕ)
                    end

                    for i in 1:Hsize-1
                        Hmat[i, i+1] = -t
                        Hmat[i+1, i] = -t
                    end

                    Hmat[1, Hsize] = -t * exp(-im * k1)
                    Hmat[Hsize, 1] = -t * exp(im * k1)

                    Hmat
                end
                @test H₀((0.0, 0.0), (6, 1)) == Flux2d((0.0, 0.0), (6, 1))
                H(k) = H₀(k, (6, 1))

                N = 51
                k = range(-π, π, length=N)
                # bandsum = (-8027.411637742611, -3926.813524046014, -1090.4729277337524, 1090.4729277337524, 3926.813524046014, 8027.411637742611)
                bandsum = (-8026.922381020279, -3931.320415546371, -1076.7367579094318, 1076.7367579094316, 3931.320415546371, 8026.92238102028)
                result = showBand(H)

                @test result.k[:, 1] == k
                @test result.k[:, 2] == k
                for i in 1:size(result.Ene, 3)
                    @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
                end

                @test abs(sum(result.Ene)) < 1e-10

                @test calcChern(H) == (TopologicalNumber=[1, 1, -2, -2, 1, 1], Total=0)

                prob = FCProblem(H)
                @test solve(prob).TopologicalNumber == [1, 1, -2, -2, 1, 1]

                C1 = zeros(6)
                C2 = zeros(6)
                N = 51
                for j in 1:N
                    for i in 1:N
                        C1 .+= calcBerryFlux(H, [i - 1, j - 1]).TopologicalNumber
                        C2 .+= calcBerryFlux(H, [i - 1, j - 1], rounds=false).TopologicalNumber
                    end
                end

                @test C1 == calcChern(H).TopologicalNumber
                @test C2 ≈ [1, 1, -2, -2, 1, 1]


                H(k, p) = H₀(k, (6, p))

                param = 1:6 # range(1, 6, length=6)
                result = calcPhaseDiagram(H, param, "Chern")
                calcPhaseDiagram(H, param, "Chern"; progress=true)
                result_MPI = calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI))
                calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI), progress=true)
                @test result.nums == result_MPI.nums

                @test result.nums[1, :] == [1, 1, -2, -2, 1, 1]
                @test result.nums[2, :] == [0, 1, -1, -1, 1, 0]
                @test result.nums[3, :] == [0, 0, 0, 0, 0, 0]
                @test result.nums[4, :] == [0, -1, 1, 1, -1, 0]
                @test result.nums[5, :] == [-1, -1, 2, 2, -1, -1]
                @test result.nums[6, :] == [0, 0, 0, 0, 0, 0]
                # @test result.nums == num

                prob = FCProblem(H)
                result = calcPhaseDiagram(prob, param)
                @test result.nums[1, :] == [1, 1, -2, -2, 1, 1]
                @test result.nums[2, :] == [0, 1, -1, -1, 1, 0]
                @test result.nums[3, :] == [0, 0, 0, 0, 0, 0]
                @test result.nums[4, :] == [0, -1, 1, 1, -1, 0]
                @test result.nums[5, :] == [-1, -1, 2, 2, -1, -1]
                @test result.nums[6, :] == [0, 0, 0, 0, 0, 0]

                fig = plot1D(result.nums, result.param; disp=false)
                @test typeof(fig) == Figure
                plotclose()

                fig = plot1D(result; disp=false)
                @test typeof(fig) == Figure
                plotclose()

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
                @test H₀((0.0, 0.0), (0.5, 1.0)) == Haldane((0.0, 0.0), (0.5, 1.0))
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

                prob = FCProblem(H)
                @test solve(prob).TopologicalNumber == [-1, 1]


                H(k, p) = H₀(k, (p, 2.5))

                param = range(-π, π, length=10)
                result = calcPhaseDiagram(H, param, "Chern")
                calcPhaseDiagram(H, param, "Chern"; progress=true)
                result_MPI = calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI))
                calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI), progress=true)
                @test result.nums == result_MPI.nums

                num = [0 0; 1 -1; 1 -1; 1 -1; 0 0; 0 0; -1 1; -1 1; -1 1; 0 0]
                @test result.nums == num

                prob = FCProblem(H)
                @test calcPhaseDiagram(prob, param).nums == num

                fig = plot1D(result.nums, result.param; disp=false)
                @test typeof(fig) == Figure
                plotclose()

                fig = plot1D(result; disp=false)
                @test typeof(fig) == Figure
                plotclose()


                result = calcPhaseDiagram(H, param, "Chern"; rounds=false)
                num = [1.1109801146677412e-15 -2.850714348572913e-16; 0.9999999999999996 -0.9999999999999991; 1.0000000000000007 -1.0000000000000002; 1.0 -1.000000000000001; 2.196996748971226e-15 4.455768167494754e-18; 6.938689780587844e-16 1.9512837637595136e-16; -0.999999999999998 0.9999999999999998; -0.9999999999999998 1.0; -1.0 0.9999999999999997; 6.668909048176786e-16 -2.7518073289070415e-16]
                @test result.nums ≈ num

                # param1 = range(-π, π, length=6)
                param1 = range(-π, π, length=7)
                # param2 = range(-6.0, 6.0, length=6)
                param2 = [-6.0, -3.6, -1.2, 1.2, 3.6, 6.0]
                result = calcPhaseDiagram(H₀, param1, param2, "Chern")
                num = zeros(2, 7, 6)
                # num[:, :, 1] = [0 0 0 0 0 0; 0 0 0 0 0 0]
                # num[:, :, 2] = [0 1 0 0 -1 0; 0 -1 0 0 1 0]
                # num[:, :, 3] = [0 1 1 -1 -1 0; 0 -1 -1 1 1 0]
                # num[:, :, 4] = [0 1 1 -1 -1 0; 0 -1 -1 1 1 0]
                # num[:, :, 5] = [0 1 0 0 -1 0; 0 -1 0 0 1 0]
                # num[:, :, 6] = [0 0 0 0 0 0; 0 0 0 0 0 0]
                num[:, :, 1] = [0 0 0 0 0 0 0; 0 0 0 0 0 0 0]
                num[:, :, 2] = [0 1 1 0 -1 -1 0; 0 -1 -1 0 1 1 0]
                num[:, :, 3] = [0 1 1 0 -1 -1 0; 0 -1 -1 0 1 1 0]
                num[:, :, 4] = [0 1 1 0 -1 -1 0; 0 -1 -1 0 1 1 0]
                num[:, :, 5] = [0 1 1 0 -1 -1 0; 0 -1 -1 0 1 1 0]
                num[:, :, 6] = [0 0 0 0 0 0 0; 0 0 0 0 0 0 0]
                # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
                @test result.nums == num

                prob = FCProblem(H₀)
                @test calcPhaseDiagram(prob, param1, param2).nums == num

                fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
                @test typeof(fig) == Figure
                plotclose()

                fig = plot2D(result; disp=false)
                @test typeof(fig) == Figure
                plotclose()

                param1 = range(-π, π, length=6)
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
                t, λₛₒ = p

                R1 = 0
                R2 = 0
                R3 = 2λₛₒ * (sin(k1) - sin(k2) - sin(k1 - k2))
                R4 = -t * (sin(k1) + sin(k2))
                R0 = -t * (cos(k1) + cos(k2) + 1)

                s0 = [1 0; 0 1]
                sx = [0 1; 1 0]
                sy = [0 -im; im 0]
                sz = [1 0; 0 -1]

                a1 = kron(sz, sx)
                a2 = kron(sz, sy)
                a3 = kron(sz, sz)
                a4 = kron(sy, s0)
                a0 = kron(sx, s0)

                R1 .* a1 .+ R2 .* a2 .+ R3 .* a3 .+ R4 .* a4 .+ R0 .* a0
            end
            @test H₀((0.0, 0.0), (0.5, 1.0)) == KaneMele((0.0, 0.0), (0.5, 1.0))
            H(k) = H₀(k, (1.0, 1.0))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-7404.378662190171, -7404.378662190167, 7404.378662190169, 7404.37866219017)
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


            prob = Z2Problem(H)
            @test solve(prob).TopologicalNumber == [1, 1]


            H(k, p) = H₀(k, (p, 1.0))

            param = range(-2.0, 2.0, length=11)
            result = calcPhaseDiagram(H, param, "Z2")
            calcPhaseDiagram(H, param, "Z2"; progress=true)
            result_MPI = calcPhaseDiagram(H, param, "Z2"; parallel=UseMPI(MPI))
            calcPhaseDiagram(H, param, "Z2"; parallel=UseMPI(MPI), progress=true)
            @test result.nums == result_MPI.nums

            num = [1 1; 1 1; 1 1; 1 1; 1 1; 0 0; 1 1; 1 1; 1 1; 1 1; 1 1]
            @test result.nums == num

            prob = Z2Problem(H)
            @test calcPhaseDiagram(prob, param).nums == num

            fig = plot1D(result.nums, result.param; disp=false)
            @test typeof(fig) == Figure
            plotclose()

            fig = plot1D(result; disp=false)
            @test typeof(fig) == Figure
            plotclose()

            result = calcPhaseDiagram(H, param, "Z2"; rounds=false)
            num = [1.000000000000001 0.9999999999999991; 0.9999999999999996 1.0; 1.0 1.0000000000000009; 0.9999999999999996 1.0; 0.9999999999999996 1.0; 0.0 0.0; 0.9999999999999996 1.0; 0.9999999999999996 1.0; 1.0 1.0000000000000009; 0.9999999999999996 1.0; 1.000000000000001 0.9999999999999991]
            @test result.nums ≈ num


            param = range(-2.0, 2.0, length=3)
            result = calcPhaseDiagram(H₀, param, param, "Z2")
            num = zeros(2, 3, 3)
            num[:, :, 1] = [1 0 1; 1 0 1]
            num[:, :, 2] = [1 0 1; 1 0 1]
            num[:, :, 3] = [1 0 1; 1 0 1]
            # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
            @test result.nums == num

            prob = Z2Problem(H₀)
            @test calcPhaseDiagram(prob, param, param).nums == num

            fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
            @test typeof(fig) == Figure
            plotclose()

            fig = plot2D(result; disp=false)
            @test typeof(fig) == Figure
            plotclose()

            result = calcPhaseDiagram(H₀, param, param, "Z2"; rounds=false)
            num = zeros(2, 3, 3)
            num[:, :, 1] = [1.0 0.0 1.0; 1.0 0.0 1.0]
            num[:, :, 2] = [1.0 0.0 1.0; 1.0 0.0 1.0]
            num[:, :, 3] = [1.0 0.0 1.0; 1.0 0.0 1.0]
            @test result.nums ≈ num

        end
    end

    @testset "3D case" begin
        function H₀(k) # Weyl
            k1, k2, k3 = k
            t1 = 1
            t2 = 1
            t3 = 1
            m = 2

            h0 = 0
            hx = 2t1 * (cos(k1) - cos(2pi * 2 / 5)) + m * (2 - cos(k2 + 2pi * 1e-3) - cos(k3 + 2pi * 1e-3))
            hy = 2t2 * sin(k2 + 2pi * 1e-3)
            hz = 2t3 * sin(k3 + 2pi * 1e-3)

            s0 = [1 0; 0 1]
            sx = [0 1; 1 0]
            sy = [0 -im; im 0]
            sz = [1 0; 0 -1]

            h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
        end

        @test calcWeylNode(H₀, [4, 10, 10]; N=11) == (TopologicalNumber=[1, -1], n=[4, 10, 10], N=11)

        N = 11
        nodes = zeros(N, N, N, 2)
        for i in 1:N, j in 1:N, k in 1:N
            nodes[i, j, k, :] = calcWeylNode(H₀, [i - 1, j - 1, k - 1]; N=N, rounds=false).TopologicalNumber
        end
        Chern_i = [[round(Int, sum(nodes[i, :, :, 1])) for i in 1:N] [round(Int, sum(nodes[i, :, :, 2])) for i in 1:N]]
        @test Chern_i[:, 1] == [0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0]
        @test Chern_i[:, 1] == -Chern_i[:, 2]

        result = calcChernSurface(H₀, "k1"; kn_mesh=11)
        @test result.kn == "k1"
        @test result.nums[:, 1] == [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
        @test result.nums[:, 2] == -result.nums[:, 1]

        result = calcChernSurface(H₀, "k2"; kn_mesh=11)
        @test result.kn == "k2"
        @test result.nums[:, 1] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        @test result.nums[:, 2] == -result.nums[:, 1]

        result = calcChernSurface(H₀, "k3"; kn_mesh=11)
        @test result.kn == "k3"
        @test result.nums[:, 1] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        @test result.nums[:, 2] == -result.nums[:, 1]

        result = findWeylPoint(H₀)
        @test result.WeylPoint == [[[4000, 9990, 9990], [6000, 9990, 9990]], [[4000, 9990, 9990], [6000, 9990, 9990]]]
        @test result.Nodes == [[1, -1], [-1, 1]]
    end

    @testset "4D case" begin
        @testset "SecondChern" begin
            @testset "Lattice Dirac model" begin

                function H₀(k, p) # landau
                    k1, k2, k3, k4 = k
                    m = p

                    # Define Pauli matrices and Gamma matrices
                    σ₀ = [1 0; 0 1]
                    σ₁ = [0 1; 1 0]
                    σ₂ = [0 -im; im 0]
                    σ₃ = [1 0; 0 -1]
                    g1 = kron(σ₁, σ₀)
                    g2 = kron(σ₂, σ₀)
                    g3 = kron(σ₃, σ₁)
                    g4 = kron(σ₃, σ₂)
                    g5 = kron(σ₃, σ₃)

                    h1 = m + cos(k1) + cos(k2) + cos(k3) + cos(k4)
                    h2 = sin(k1)
                    h3 = sin(k2)
                    h4 = sin(k3)
                    h5 = sin(k4)

                    # Update the Hamiltonian matrix in place
                    h1 .* g1 .+ h2 .* g2 .+ h3 .* g3 .+ h4 .* g4 .+ h5 .* g5
                end
                @test H₀((0.0, 0.0, 0.0, 0.0), -3.0) == LatticeDirac((0.0, 0.0, 0.0, 0.0), -3.0)
                H(k) = H₀(k, -3.0)

                N = (10, 10, 10, 10)

                @test calcSecondChern(H; N).TopologicalNumber ≈ 0.8309301430562057
                @test calcSecondChern(H; N, parallel=UseMPI(MPI)).TopologicalNumber ≈ 0.8309301430562057

                prob = SCProblem(H, N)
                @test solve(prob).TopologicalNumber ≈ 0.8309301430562057

                param = range(-4.9, 4.9, length=4)
                result = calcPhaseDiagram(H₀, param, FHS2(); N=10)
                calcPhaseDiagram(H₀, param, FHS2(); N=10, progress=true)
                result_MPI = calcPhaseDiagram(H₀, param, FHS2(); N=10, parallel=UseMPI(MPI))
                calcPhaseDiagram(H₀, param, FHS2(); N=10, parallel=UseMPI(MPI), progress=true)
                @test result.nums ≈ result_MPI.nums

                nums = [0.0010237313095167225, -2.0667333080974735, 2.1572606447321454, -0.0009805850180973213]
                @test result.nums ≈ nums

                prob = SCProblem(H₀, N)
                @test calcPhaseDiagram(prob, param).nums ≈ nums


                result = calcPhaseDiagram(H₀, param, FHS2(); N=10, rounds=false)
                nums = ComplexF64[0.0010237313095167225+8.29577265997263e-17im, -2.0667333080974735-1.6878307095102013e-16im, 2.1572606447321454+3.05582947604634e-16im, -0.0009805850180973213+3.3561217905699694e-17im]
                @test result.nums ≈ nums
            end

        end
    end

    @testset "model" begin
        @testset "SSH" begin
            H(k) = SSH(k, (0.3, 0.5))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-27.505964588866973, 27.505964588866973)
            result = showBand(H)

            @test result.k == k
            for i in 1:size(result.Ene, 2)
                @test sum(result.Ene[:, i]) ≈ bandsum[i]
            end

            @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "KitaevChain" begin
            H(k) = KitaevChain(k, (0.5, 0.2))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-70.31377676120336, 70.31377676120336)
            result = showBand(H)

            @test result.k == k
            for i in 1:size(result.Ene, 2)
                @test sum(result.Ene[:, i]) ≈ bandsum[i]
            end

            @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "Flux2d" begin
            H(k) = Flux2d(k, (6, 1))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-8026.922381020279, -3931.320415546371, -1076.736757909432, 1076.7367579094318, 3931.320415546371, 8026.92238102028)
            result = showBand(H)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end
            @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "Haldane" begin
            H(k) = Haldane(k, (0.5, 0.5))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-5563.582717519327, 5209.0393625156175)
            result = showBand(H)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end
            # @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "KitaevHoneycomb" begin
            H(k) = KitaevHoneycomb(k, (1.0, 0.2))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-4136.787430207456, 4136.787430207456)
            result = showBand(H)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end
            @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "KaneMele" begin
            H(k) = KaneMele(k, (1.0, 0.3))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-4652.438396310186, -4652.438396310185, 4652.438396310186, 4652.438396310186)
            result = showBand(H)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end
            @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "BHZ" begin
            H(k) = BHZ(k, (0.5, 1.0))

            N = 51
            k = range(-π, π, length=N)
            bandsum = (-2874.262392663075, -2874.2623926630727, 8484.262392663073, 8484.262392663073)
            result = showBand(H)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end
            # @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "LatticeDirac" begin
            H(k) = LatticeDirac(k, -3.0)

            N = 11
            k = range(-π, π, length=N)
            bandsum = (-54020.051075291514, -54020.0510752915, 54020.0510752915, 54020.051075291514)
            result = showBand(H; N=N)

            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            @test result.k[:, 3] == k
            @test result.k[:, 4] == k
            for i in 1:size(result.Ene, 5)
                @test sum(result.Ene[:, :, :, :, i]) ≈ bandsum[i]
            end

            @test abs(sum(result.Ene)) < 1e-10
        end
    end

end

