abstract type TopologicalNumbersParallel end
abstract type TopologicalNumbersSingleProcess <: TopologicalNumbersParallel end

struct UseSingleThread <: TopologicalNumbersSingleProcess end


abstract type TopologicalNumbersMultiProcess <: TopologicalNumbersParallel end

struct UseThreads <: TopologicalNumbersMultiProcess end
struct UseFLoop <: TopologicalNumbersMultiProcess end
struct UseMPI <: TopologicalNumbersMultiProcess end