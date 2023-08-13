@doc raw"""
    plot1D(nums::Matrix, param_range::T; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T<:AbstractVector}
"""
function plot1D(nums::Matrix, param_range::T; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T<:AbstractVector}

    fig = Figure()

    if labels == true
        ax = Axis(fig[1, 1], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="p", ylabel="ν")
    else
        ax = Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
    end


    for i in 1:size(nums, 2)
        lines!(ax, param_range, nums[:, i], label="Band$(i)", linewidth=5)
    end
    Legend(fig[1, 2], ax)

    p = (; png, pdf, svg, filename)
    output(fig, p)
    display(fig)
    fig
end

@doc raw"""
    plot2D(nums::T1, param_range1::T2, param_range2::T2; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractArray,T2<:AbstractVector}
"""
function plot2D(nums::T1, param_range1::T2, param_range2::T2; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractArray,T2<:AbstractVector}

    fig = Figure()

    if labels == true
        ax = Axis(fig[1, 1], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(2), xlabel="p₁", ylabel="p₂")
    else
        ax = Axis(fig[1, 1], xlabelvisible=false, ylabelvisible=false)
    end


    hm = heatmap!(ax, param_range1, param_range2, nums, colormap=:jet1)
    ax.aspect = AxisAspect(1)
    Colorbar(fig[1, 2], hm)

    p = (; png, pdf, svg, filename)
    output(fig, p)
    display(fig)
    fig
end

function output(fig, p)
    @unpack png, pdf, svg, filename = p
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
end