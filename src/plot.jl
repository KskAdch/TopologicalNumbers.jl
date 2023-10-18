@doc raw"""
    plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}
"""
function plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}

    fig = Figure()

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

    fig = Figure()

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

    fig = Figure()

    if labels == true
        ax = Axis(fig[1, 1], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="p₁", ylabel="p₂")
    else
        ax = Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
    end


    hm = heatmap!(ax, param_range1, param_range2, nums, colormap=:jet1)
    ax.aspect = AxisAspect(1)
    Colorbar(fig[1, 2], hm)

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    # display(fig)
    fig
end

@doc raw"""
    plot2D(result::NamedTuple; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")
"""
function plot2D(result::NamedTuple; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")

    fig = Figure()

    if labels == true
        ax = Axis(fig[1, 1], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="p₁", ylabel="p₂")
    else
        ax = Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
    end

    nums_half = sum(@view(result.nums[1:end÷2, :, :]), dims=1)[1, :, :]

    hm = heatmap!(ax, result.param1, result.param2, nums_half, colormap=:jet1)
    ax.aspect = AxisAspect(1)
    Colorbar(fig[1, 2], hm)

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    # display(fig)
    fig
end

function output(fig, p)
    @unpack disp, png, pdf, svg, filename = p
    if png == true
        CairoMakie.activate!()
        save(filename * ".png", fig)
    end
    if pdf == true
        CairoMakie.activate!()
        save(filename * ".pdf", fig)
    end
    if svg == true
        CairoMakie.activate!()
        save(filename * ".svg", fig)
    end
    if disp == true
        GLMakie.activate!()
        display(fig)
    end
end