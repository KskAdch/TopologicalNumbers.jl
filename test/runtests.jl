using TopologicalNumbers
using LinearAlgebra
using Test

using MPI

using CondaPkg
# CondaPkg.add("numpy")
CondaPkg.add("pfapack")
# CondaPkg.add("scipy")


using PythonCall
using PythonPlot

# pyimport("sys").path.append(@__DIR__)
# const pf = pyimport("pfaffian")

const pf = pyimport("pfapack.pfaffian")
# const cpf = pyimport("pfapack.ctypes").pfaffian
const np = pyimport("numpy")

# using Aqua
# Aqua.test_all(TopologicalNumbers; ambiguities=false)

@testset "TopologicalNumbers.jl" begin

    # test pfaffian
    include("test_pfaffian.jl")
    @testset "1D case" begin
        function H₀(k, p)

            [
                0 p[1]+p[2]*exp(-im * k)
                p[1]+p[2]*exp(im * k) 0
            ]
        end
        # @test H₀(0.0, (1.0, 1.0)) == SSH(0.0, (1.0, 1.0))
        @test H₀(0.0, (1.0, 1.0)) == SSH(0.0, 1.0)
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

        BPProblem(H, 41)
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

                FCProblem(H, 41)
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

                @test C1 == [1, 1, -2, -2, 1, 1]
                @test C2 ≈ [1, 1, -2, -2, 1, 1]


                LBFProblem(H, zeros(3), 41)
                C1 = zeros(6)
                C2 = zeros(6)
                N = 51
                for j in 1:N
                    for i in 1:N
                        prob = LBFProblem(H, [i - 1, j - 1])
                        C1 .+= solve(prob).TopologicalNumber
                        prob = LBFProblem(; H, n=[i - 1, j - 1], rounds=false)
                        C2 .+= solve(prob).TopologicalNumber
                    end
                end

                @test C1 == [1, 1, -2, -2, 1, 1]
                @test C2 ≈ [1, 1, -2, -2, 1, 1]


                H(k, p) = H₀(k, (6, p))

                param = 1:6 # range(1, 6, length=6)
                result = calcPhaseDiagram(H, param, "Chern")
                calcPhaseDiagram(H, param, "Chern"; progress=true, rounds = false)
                result_MPI = calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI))
                calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI), progress=true, rounds = false)
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
                    hx = J * (1 + cos(k1) + cos(k2))
                    hy = J * (-sin(k1) + sin(k2))
                    hz = M - 2K * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))

                    s0 = [1 0; 0 1]
                    sx = [0 1; 1 0]
                    sy = [0 -im; im 0]
                    sz = [1 0; 0 -1]

                    h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz
                end
                # @test H₀((0.0, 0.0), (0.5, 1.0)) == Haldane((0.0, 0.0), (0.5, 1.0))
                @test H₀((0.0, 0.0), (0.5, 1.0)) == Haldane((0.0, 0.0), (1, 0.5, 1.0))
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

                @test calcChern(H) == (TopologicalNumber=[1, -1], Total=0)

                prob = FCProblem(H)
                @test solve(prob).TopologicalNumber == [1, -1]


                H(k, p) = H₀(k, (p, 2.5))

                param = range(-π, π, length=10)
                result = calcPhaseDiagram(H, param, "Chern")
                calcPhaseDiagram(H, param, "Chern"; progress=true)
                result_MPI = calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI))
                calcPhaseDiagram(H, param, "Chern"; parallel=UseMPI(MPI), progress=true)
                @test result.nums == result_MPI.nums

                num = [0 0; -1 1; -1 1; -1 1; 0 0; 0 0; 1 -1; 1 -1; 1 -1; 0 0]
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
                num = [1.664954913440598e-15 3.642517302686141e-17; -0.9999999999999998 1.0000000000000002; -0.9999999999999998 0.9999999999999996; -0.9999999999999987 1.0000000000000002; -6.659742003501615e-16 -1.0767325801951022e-15; 6.670753241870261e-16 -7.62595798192664e-16; 1.0 -0.9999999999999994; 1.0 -0.9999999999999993; 1.0000000000000036 -0.9999999999999999; 7.767764937404728e-16 -3.7599154015552054e-18]
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
                num[:, :, 2] = [0 -1 -1 0 1 1 0; 0 1 1 0 -1 -1 0]
                num[:, :, 3] = [0 -1 -1 0 1 1 0; 0 1 1 0 -1 -1 0]
                num[:, :, 4] = [0 -1 -1 0 1 1 0; 0 1 1 0 -1 -1 0]
                num[:, :, 5] = [0 -1 -1 0 1 1 0; 0 1 1 0 -1 -1 0]
                num[:, :, 6] = [0 0 0 0 0 0 0; 0 0 0 0 0 0 0]
                # num = [0 1 0; 0 1 1;;; 0 0 0; 0 0 0;;; 0 1 0; 1 1 0]
                @test result.nums == num

                Hal(k, p) = Haldane(k, (1, p[1], p[2]))
                calcPhaseDiagram(Hal, [0.0, π], [0.0, 1.0], "Chern"; parallel=UseMPI(MPI), progress=true)

                prob = FCProblem(H₀)
                @test calcPhaseDiagram(prob, param1, param2).nums == num

                fig = plot2D(result.nums[1, :, :], result.param1, result.param2; disp=false)
                @test typeof(fig) == Figure
                plotclose()

                fig = plot2D(result; disp=false)
                @test typeof(fig) == Figure
                plotclose()

                param1 = range(-π, π, length=7)
                result = calcPhaseDiagram(H₀, param1, param2, "Chern"; rounds=false)
                num = zeros(2, 7, 6)
                num[:, :, 1] = [-1.0627809926925885e-15 -3.0679105522077703e-16 -3.156342701459608e-16 -1.0811601221010123e-15 -3.0630256200953273e-16 -2.0039600971368465e-16 -1.057245715422552e-15; -2.2168547246893346e-16 2.2205431120762855e-16 2.2195077752659136e-16 -2.214697773001059e-16 4.4363085761630395e-16 4.437834619482495e-16 -2.2168547246893346e-16]
                num[:, :, 2] = [-1.051257930787215e-15 -0.9999999999999991 -0.9999999999999996 -1.0069835284171914e-15 1.0000000000000004 1.0 -1.0402089457640254e-15; 1.554088989975483e-15 0.999999999999999 0.9999999999999983 1.3327389234940762e-15 -0.9999999999999987 -0.9999999999999987 1.554088989975483e-15]
                num[:, :, 3] = [1.8983141818413044e-16 -0.9999999999999997 -1.0000000000000004 6.470949266457012e-16 0.9999999999999996 0.9999999999999994 1.8958121178829055e-16; 1.2191110233921236e-15 0.9999999999999999 0.9999999999999992 1.2217619170170142e-15 -0.9999999999999987 -0.9999999999999998 9.970664184670923e-16]
                num[:, :, 4] = [9.997173120919829e-16 -0.9999999999999998 -0.9999999999999991 5.529772086170295e-16 1.0000000000000002 0.9999999999999999 1.2217619170170142e-15; 1.4124356215355783e-16 0.9999999999999994 0.9999999999999997 2.1216449596453474e-16 -1.0000000000000004 -1.0000000000000007 1.556563133346148e-16]
                num[:, :, 5] = [1.3327389234940762e-15 -0.9999999999999978 -0.9999999999999978 1.554088989975483e-15 0.9999999999999981 0.999999999999999 1.3327389234940762e-15; -1.0139483254186342e-15 1.0 1.0000000000000013 -1.0376303100206914e-15 -0.9999999999999996 -0.9999999999999996 -1.0088083095454722e-15]
                num[:, :, 6] = [-2.214697773001059e-16 4.4363085761630395e-16 4.437834619482495e-16 -2.2168547246893346e-16 2.2205431120762855e-16 2.2195077752659136e-16 -2.214697773001059e-16; -1.0723230910341474e-15 -1.858487186710715e-16 -3.3489079968593663e-16 -1.0583506139248706e-15 -3.0624356021802653e-16 -3.1667152799143036e-16 -1.0811612005768565e-15]
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
            # @test H₀((0.0, 0.0), (0.5, 1.0)) == KaneMele((0.0, 0.0), (0.5, 1.0))
            @test H₀((0.0, 0.0), (1, 1.0)) == KaneMele((0.0, 0.0), 1.0)
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

            @test calcZ2(H) == (TopologicalNumber=[1, 1], TRTopologicalNumber=nothing, Total=0)

            @test norm(calcZ2(H, rounds=false).TopologicalNumber - calcZ2(H).TopologicalNumber) < 1e-10

            @test calcZ2(H, TR=true) == (TopologicalNumber=[1, 1], TRTopologicalNumber=[1, 1], Total=0)

            @test norm(calcZ2(H, rounds=false, TR=true).TRTopologicalNumber - calcZ2(H, TR=true).TRTopologicalNumber) < 1e-10


            Z2Problem(H, 2)
            Z2Problem(H, 2, 40)
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

            param1 = range(-1.0, 1.0, length = 2)
            param2 = range(-0.5, 0.5, length = 2)
            res = calcPhaseDiagram(BHZ, param1, param2, "Z2"; parallel=UseMPI(MPI), progress=true)
            # res.nums .≈ [0.9999999999999989 1.0; 0.9999999999999987 0.9999999999999993;;; 1.0 0.9999999999999989; 0.9999999999999991 0.9999999999999987]
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

        WNProblem(H₀, [4, 10, 10])
        prob = WNProblem(H₀, [4, 10, 10], 11)
        result = solve(prob)
        @test result.TopologicalNumber == [1, -1]
        @test result.n == [4, 10, 10]
        @test result.N == 11

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


        WCSProblem(H₀, "k1")
        WCSProblem(H₀, "k1", 11, 41)
        prob = WCSProblem(H₀, "k1", 11)
        result = solve(prob)
        @test result.kn == "k1"
        @test result.nums[:, 1] == [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
        @test result.nums[:, 2] == -result.nums[:, 1]

        prob = WCSProblem(H₀, "k2", 11)
        result = solve(prob)
        @test result.kn == "k2"
        @test result.nums[:, 1] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        @test result.nums[:, 2] == -result.nums[:, 1]

        prob = WCSProblem(H₀, "k3", 11)
        result = solve(prob)
        @test result.nums[:, 1] == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        @test result.nums[:, 2] == -result.nums[:, 1]


        result = findWeylPoint(H₀)
        @test result.WeylPoint == [[[4000, 9990, 9990], [6000, 9990, 9990]], [[4000, 9990, 9990], [6000, 9990, 9990]]]
        @test result.Nodes == [[1, -1], [-1, 1]]

        WPProblem(H₀, 41)
        prob = WPProblem(H₀)
        result = solve(prob)
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

                SCProblem(H)
                SCProblem(H, 1)
                prob = SCProblem(; H, N)
                @test solve(prob).TopologicalNumber ≈ 0.8309301430562057

                H(k, p) = LatticeDirac(k, p[1]) + p[2]*Matrix{Float64}(I, 4, 4)
                prob = SCProblem(H, 1, 31)
                param1 = [-1.0, 1.0]
                param2 = [-0.5, 0.5]
                res = calcPhaseDiagram(prob, param1, param2)
                @test round.(res.nums, digits = 20) ≈ round.([2.3975776848447665e-29 2.57553745875129e-29; 2.3675737006667523e-29 2.5844675711354074e-29], digits = 20)

                param = range(-4.9, 4.9, length=4)
                result = calcPhaseDiagram(H₀, param, FHS2(); N=10)
                calcPhaseDiagram(H₀, param, FHS2(); N=10, progress=true)
                result_MPI = calcPhaseDiagram(H₀, param, FHS2(); N=10, parallel=UseMPI(MPI))
                calcPhaseDiagram(H₀, param, FHS2(); N=10, parallel=UseMPI(MPI), progress=true)
                @test result.nums ≈ result_MPI.nums

                nums = [0.0010237313095167225, -2.0667333080974735, 2.1572606447321454, -0.0009805850180973213]
                @test result.nums ≈ nums

                prob = SCProblem(; H=H₀, N)
                @test calcPhaseDiagram(prob, param).nums ≈ nums


                result = calcPhaseDiagram(H₀, param, FHS2(); N=10, rounds=false)
                nums = ComplexF64[0.0010237313095167225+8.29577265997263e-17im, -2.0667333080974735-1.6878307095102013e-16im, 2.1572606447321454+3.05582947604634e-16im, -0.0009805850180973213+3.3561217905699694e-17im]
                @test result.nums ≈ nums

                fig = plot1D(FHS2(), result.nums, result.param; disp=false)
                @test typeof(fig) == Figure
                plotclose()
            end

        end
    end

    @testset "model" begin
        @testset "SSH" begin
            # H(k) = SSH(k, (0.3, 0.5))
            H(k) = SSH(k, 0.6)

            N = 51
            k = range(-π, π, length=N)
            # bandsum = (-27.505964588866973, 27.505964588866973)
            bandsum = (-55.01192917773395, 55.01192917773395)
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
            # H(k) = Haldane(k, (0.5, 0.5))
            H(k) = Haldane(k, (1, 0.5, 0.5))

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
            # H(k) = KitaevHoneycomb(k, (1.0, 0.2))
            H(k) = KitaevHoneycomb(k, 0.2)

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

        @testset "ThoulessPump" begin
            H(k) = ThoulessPump(k, (-1.0, 0.5))
        
            N = 51
            k = range(-π, π, length=N)
            bandsum = (-6140.668216638056, -3853.1610208575285, 3853.1610208575285, 6140.668216638056)
            result = showBand(H)
        
            @test result.k[:, 1] == k
            @test result.k[:, 2] == k
            for i in 1:size(result.Ene, 3)
                @test sum(result.Ene[:, :, i]) ≈ bandsum[i]
            end
            @test abs(sum(result.Ene)) < 1e-10
        end

        @testset "KaneMele" begin
            # H(k) = KaneMele(k, (1.0, 0.3))
            H(k) = KaneMele(k, 0.3)

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
            # bandsum = (-2874.262392663075, -2874.2623926630727, 8484.262392663073, 8484.262392663073)
            bandsum = (-4276.762392663073, -4276.762392663073, 7081.762392663072, 7081.762392663074)
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

