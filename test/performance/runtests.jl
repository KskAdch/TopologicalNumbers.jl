using Test
using TopologicalNumbers

const K2 = (0.37, 1.11)
const K4 = (0.17, 0.43, 0.79, 1.23)
const HALDANE_PARAMETERS = (1.0, 0.5, 0.4)

# コンパイル時間を除外し、実行時間とヒープ割り当て量を同じ条件で測定する。
function measure_performance(f)
    result = f()
    GC.gc()
    allocated = @allocated result = f()
    GC.gc()
    elapsed = @elapsed result = f()
    return (; result, allocated, elapsed)
end

format_mib(bytes) = round(bytes / 2^20; digits=3)
format_ms(seconds) = round(seconds * 1_000; digits=3)

@testset "Performance regression" begin
    @testset "Preset Hamiltonians" begin
        haldane = measure_performance(() -> Haldane(K2, HALDANE_PARAMETERS))
        dirac = measure_performance(() -> LatticeDirac(K4, -3.0))

        @info "Haldane" elapsed_ms = format_ms(haldane.elapsed) allocated_mib = format_mib(
            haldane.allocated
        )
        @info "LatticeDirac" elapsed_ms = format_ms(dirac.elapsed) allocated_mib = format_mib(
            dirac.allocated
        )

        @test size(haldane.result) == (2, 2)
        @test size(dirac.result) == (4, 4)
        @test haldane.allocated <= 16 * 2^10
        @test dirac.allocated <= 64 * 2^10
    end

    @testset "First Chern number" begin
        hamiltonian(k) = Haldane(k, HALDANE_PARAMETERS)
        chern = measure_performance(() -> calcChern(hamiltonian; N=11))

        @info "calcChern (N=11)" elapsed_ms = format_ms(chern.elapsed) allocated_mib = format_mib(
            chern.allocated
        )

        @test chern.result == (TopologicalNumber=[1, -1], Total=0)
        @test chern.allocated <= 16 * 2^20
    end
end
