abstract type TopologicalNumbersParallel end
abstract type TopologicalNumbersSingleProcess <: TopologicalNumbersParallel end

@doc raw"""
    UseSingleThread <: TopologicalNumbersSingleProcess

A struct representing the use of a single thread for parallel processing in the TopologicalNumbers module.

"""
struct UseSingleThread <: TopologicalNumbersSingleProcess end

abstract type TopologicalNumbersMultiProcess <: TopologicalNumbersParallel end

# struct UseThreads <: TopologicalNumbersMultiProcess end

@doc raw"""
    UseMPI{T}(MPI::T) <: TopologicalNumbersMultiProcess

A struct representing the use of MPI for parallel computing in the TopologicalNumbers package.

# Arguments
- `MPI`: An object representing the MPI library.
"""
Base.@kwdef struct UseMPI{T} <: TopologicalNumbersMultiProcess
    MPI::T
end
