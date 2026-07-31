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

const _PYTHONPLOT_PACKAGE = Base.PkgId(
    Base.UUID("274fc56d-3b97-40fa-a1cd-1b4a50311bf9"), "PythonPlot"
)

"""
利用者が読み込んだ PythonPlot モジュールを返す。

数値計算だけを行う利用者には Python 環境が不要なため、パッケージ読込時には初期化しない。
"""
function _pythonplot()
    pythonplot = get(Base.loaded_modules, _PYTHONPLOT_PACKAGE, nothing)
    if isnothing(pythonplot)
        throw(
            ArgumentError(
                "描画機能には PythonPlot.jl が必要です。" *
                "`import Pkg; Pkg.add(\"PythonPlot\"); using PythonPlot` を実行してください。",
            ),
        )
    end
    return pythonplot
end
