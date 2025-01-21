using Accessors: @reset
using Distributed
using LaTeXStrings
using LinearAlgebra
using Parameters: @unpack
using ProgressBars
using SparseArrays
using StaticArrays

# For AutoMerge CI
if haskey(ENV, "GITHUB_ACTIONS") && haskey(ENV, "AUTOMERGE_GITHUB_TOKEN")
    ENV["MPLBACKEND"] = "Agg"
end
# ENV["MPLBACKEND"] = "Agg"
using PythonPlot