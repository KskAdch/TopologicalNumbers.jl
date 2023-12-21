@doc raw"""
    plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}
"""
function plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}

    fig = figure()

    if labels == true
        ax = Axis(fig[1, 1], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="p", ylabel="ν")
    else
        ax = Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
    end


    for i in 1:size(nums, 2)
        scatter!(ax, param_range, nums[:, i], label="Band$(i)", markersize=15)
    end
    Legend(fig[1, 2], ax)

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    # display(fig)
    fig
end

@doc raw"""
    plot1D(result::NamedTuple; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")
"""
function plot1D(result::NamedTuple; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")

    fig = figure()

    if labels == true
        ax = Axis(fig[1, 1], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="p", ylabel="ν")
    else
        ax = Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
    end


    for i in 1:size(result.nums, 2)
        scatter!(ax, result.param, result.nums[:, i], label="Band$(i)", markersize=15)
    end
    Legend(fig[1, 2], ax)

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    # display(fig)
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

    im = ax.imshow(nums, cmap="jet", interpolation="none", origin="lower", extent=(param_range1[1], param_range1[end], param_range2[1], param_range2[end]), aspect="auto")
    fig.colorbar(im, ax=ax, ticks=matplotlib.ticker.MaxNLocator(integer=true))


    # hm = heatmap!(ax, param_range1, param_range2, nums, colormap=:jet1)
    # ax.aspect = AxisAspect(1)
    # Colorbar(fig[1, 2], hm)

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    # display(fig)
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

    im = ax.imshow(nums_half, cmap="jet", interpolation="none", origin="lower", extent=(result.param1[1], result.param1[end], result.param2[1], result.param2[end]), aspect="auto")
    fig.colorbar(im, ax=ax, ticks=matplotlib.ticker.MaxNLocator(integer=true))

    # hm = heatmap!(ax, result.param1, result.param2, nums_half, colormap=:jet1)
    # ax.aspect = AxisAspect(1)
    # Colorbar(fig[1, 2], hm)

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    # display(fig)
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