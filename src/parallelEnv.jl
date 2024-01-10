abstract type TopologicalNumbersParallel end
abstract type TopologicalNumbersSingleProcess <: TopologicalNumbersParallel end

@doc raw"""
"""
struct UseSingleThread <: TopologicalNumbersSingleProcess end


abstract type TopologicalNumbersMultiProcess <: TopologicalNumbersParallel end

# struct UseThreads <: TopologicalNumbersMultiProcess end

@doc raw"""
"""
Base.@kwdef struct UseMPI{T} <: TopologicalNumbersMultiProcess
    MPI::T
end