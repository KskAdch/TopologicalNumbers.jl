abstract type TopologicalNumbersParallel end
abstract type TopologicalNumbersSingleProcess <: TopologicalNumbersParallel end

struct UseSingleThread <: TopologicalNumbersSingleProcess end


abstract type TopologicalNumbersMultiProcess <: TopologicalNumbersParallel end

struct UseThreads <: TopologicalNumbersMultiProcess end
Base.@kwdef struct UseMPI{T} <: TopologicalNumbersMultiProcess
    MPI::T
end