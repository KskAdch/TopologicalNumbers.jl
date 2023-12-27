function marker(i)
    markers = ("o", "v", "^", ",", "<", ">", "8", "s", "*", "h", "H", "D", "d", ".", "1", "2", "3", "4", "+", "x", "|", "_", "p")
    if i > length(markers)
        i = i % length(markers)
    end
    return markers[i]
end

@doc raw"""
    plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}
"""
function plot1D(::T1, nums::T2, param_range::T3; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:SecondChernAlgorithms,T2<:AbstractVector,T3<:AbstractVector}

    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p")
        ax.set_ylabel(L"\nu")
    end
    ax.grid()
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=true))

    if disp == true || png == true || pdf == true || svg == true
        ax.scatter(param_range, nums, marker=marker(i), label="Band$(i)")
        ax.legend()
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    fig
end

@doc raw"""
    plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}
"""
function plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}

    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p")
        ax.set_ylabel(L"\nu")
    end
    ax.grid()
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=true))

    if disp == true || png == true || pdf == true || svg == true
        for i in axes(nums, 2)
            ax.scatter(param_range, nums[:, i], marker=marker(i), label="Band$(i)")
        end
        ax.legend()
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    fig
end

@doc raw"""
    plot1D(result::NamedTuple; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")
"""
function plot1D(result::NamedTuple; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")

    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p")
        ax.set_ylabel(L"\nu")
    end
    ax.grid()
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=true))

    if disp == true || png == true || pdf == true || svg == true
        for i in axes(result.nums, 2)
            ax.scatter(result.param, result.nums[:, i], marker=marker(i), label="Band$(i)")
        end
        ax.legend()
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    fig
end

@doc raw"""
    plot2D(nums::T1, param_range1::T2, param_range2::T3; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractArray,T2<:AbstractVector,T3<:AbstractVector}
"""
function plot2D(nums::T1, param_range1::T2, param_range2::T3; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractArray,T2<:AbstractVector,T3<:AbstractVector}

    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p_1")
        ax.set_ylabel(L"p_2")
    end

    if disp == true || png == true || pdf == true || svg == true
        im = ax.imshow(nums, cmap="jet", interpolation="none", origin="lower", extent=(param_range1[1], param_range1[end], param_range2[1], param_range2[end]), aspect="auto")
        fig.colorbar(im, ax=ax, ticks=matplotlib.ticker.MaxNLocator(integer=true))
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    fig
end

@doc raw"""
    plot2D(result::NamedTuple; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")
"""
function plot2D(result::NamedTuple; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")

    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p_1")
        ax.set_ylabel(L"p_2")
    end

    nums_half = transpose(sum(@view(result.nums[1:end÷2, :, :]), dims=1)[1, :, :])

    if disp == true || png == true || pdf == true || svg == true
        im = ax.imshow(nums_half, cmap="jet", interpolation="none", origin="lower", extent=(result.param1[1], result.param1[end], result.param2[1], result.param2[end]), aspect="auto")
        fig.colorbar(im, ax=ax, ticks=matplotlib.ticker.MaxNLocator(integer=true))
    end


    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    fig
end

function output(fig, p)
    @unpack disp, png, pdf, svg, filename = p
    if png == true
        savefig(filename * ".png")
    end
    if pdf == true
        savefig(filename * ".pdf")
    end
    if svg == true
        savefig(filename * ".svg")
    end
    if disp == true
        plotshow()
    end
end