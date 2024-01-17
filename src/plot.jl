function marker(i)
    markers = (
        "o",
        "v",
        "^",
        ",",
        "<",
        ">",
        "8",
        "s",
        "*",
        "h",
        "H",
        "D",
        "d",
        ".",
        "1",
        "2",
        "3",
        "4",
        "+",
        "x",
        "|",
        "_",
        "p",
    )
    if i > length(markers)
        i = i % length(markers)
    end
    return markers[i]
end

@doc raw"""
    plot1D(::T1, nums::T2, param_range::T3; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:SecondChernAlgorithms,T2<:AbstractVector,T3<:AbstractVector}
    
Plot a 1D phase diagram.

# Arguments
- `::T1`: An object of type `T1` that implements the `SecondChernAlgorithms` interface.
- `nums::T2`: A vector of numbers representing the values to be plotted.
- `param_range::T3`: A vector representing the parameter range.
- `labels::Bool=true`: Whether to display labels on the plot. Default is `true`.
- `disp::Bool=true`: Whether to display the plot. Default is `true`.
- `png::Bool=false`: Whether to save the plot as a PNG file. Default is `false`.
- `pdf::Bool=false`: Whether to save the plot as a PDF file. Default is `false`.
- `svg::Bool=false`: Whether to save the plot as an SVG file. Default is `false`.
- `filename::String="phaseDiagram"`: The filename for the saved plot. Default is `"phaseDiagram"`.

# Examples
```julia
plot1D(obj, nums, param_range)
```

# Returns
- `fig`: The matplotlib figure object representing the plot.
"""
function plot1D(
    ::T1,
    nums::T2,
    param_range::T3;
    labels::Bool=true,
    disp::Bool=true,
    png::Bool=false,
    pdf::Bool=false,
    svg::Bool=false,
    filename::String="phaseDiagram",
) where {T1<:SecondChernAlgorithms,T2<:AbstractVector,T3<:AbstractVector}
    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p")
        ax.set_ylabel(L"\nu")
    end
    ax.grid()
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(; integer=true))

    if disp == true || png == true || pdf == true || svg == true
        ax.scatter(param_range, nums; marker=marker(1))
        # ax.legend()
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    return fig
end

@doc raw"""
    plot1D(nums::T1, param_range::T2; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractMatrix,T2<:AbstractVector}
    
Plot a 1D phase diagram.

# Arguments
- `nums::T1`: A matrix or vector containing the data to be plotted.
- `param_range::T2`: A vector representing the parameter range.
- `labels::Bool=true`: Whether to display axis labels. Default is `true`.
- `disp::Bool=true`: Whether to display the plot. Default is `true`.
- `png::Bool=false`: Whether to save the plot as a PNG file. Default is `false`.
- `pdf::Bool=false`: Whether to save the plot as a PDF file. Default is `false`.
- `svg::Bool=false`: Whether to save the plot as an SVG file. Default is `false`.
- `filename::String="phaseDiagram"`: The filename for the saved plot. Default is "phaseDiagram".

# Examples
```julia
nums = [1 2 3; 4 5 6]
param_range = [0.1, 0.2, 0.3]
plot1D(nums, param_range)
```

# Returns
- `fig`: The matplotlib figure object representing the plot.
"""
function plot1D(
    nums::T1,
    param_range::T2;
    labels::Bool=true,
    disp::Bool=true,
    png::Bool=false,
    pdf::Bool=false,
    svg::Bool=false,
    filename::String="phaseDiagram",
) where {T1<:Union{AbstractVector,AbstractMatrix},T2<:AbstractVector}
    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p")
        ax.set_ylabel(L"\nu")
    end
    ax.grid()
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(; integer=true))

    if disp == true || png == true || pdf == true || svg == true
        if nums isa AbstractVector
            ax.scatter(param_range, nums; marker=marker(1))
        else
            for i in axes(nums, 2)
                ax.scatter(param_range, nums[:, i]; marker=marker(i), label="Band$(i)")
            end
        end
        ax.legend()
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    return fig
end

@doc raw"""
    plot1D(result::NamedTuple; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")
    
Plot a 1D phase diagram.

# Arguments
- `result::NamedTuple`: A named tuple containing the result data.
- `labels::Bool=true`: Whether to display axis labels.
- `disp::Bool=true`: Whether to display the plot.
- `png::Bool=false`: Whether to save the plot as a PNG file.
- `pdf::Bool=false`: Whether to save the plot as a PDF file.
- `svg::Bool=false`: Whether to save the plot as an SVG file.
- `filename::String="phaseDiagram"`: The filename for the saved plot.

# Returns
- `fig`: The matplotlib figure object.

# Example
```julia
julia>
```
"""
function plot1D(
    result::NamedTuple;
    labels::Bool=true,
    disp::Bool=true,
    png::Bool=false,
    pdf::Bool=false,
    svg::Bool=false,
    filename::String="phaseDiagram",
)
    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p")
        ax.set_ylabel(L"\nu")
    end
    ax.grid()
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(; integer=true))

    if disp == true || png == true || pdf == true || svg == true
        if result.nums isa AbstractVector
            ax.scatter(result.param, result.nums; marker=marker(1))
        else
            for i in axes(result.nums, 2)
                ax.scatter(
                    result.param, result.nums[:, i]; marker=marker(i), label="Band$(i)"
                )
            end
        end
        ax.legend()
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    return fig
end

@doc raw"""
    plot2D(nums::T1, param_range1::T2, param_range2::T3; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram") where {T1<:AbstractArray,T2<:AbstractVector,T3<:AbstractVector}

Plot a 2D phase diagram.

# Arguments
- `nums::T1`: An array of numbers representing the phase diagram.
- `param_range1::T2`: A vector representing the range of the first parameter.
- `param_range2::T3`: A vector representing the range of the second parameter.
- `labels::Bool=true`: Whether to display axis labels. Default is `true`.
- `disp::Bool=true`: Whether to display the plot. Default is `true`.
- `png::Bool=false`: Whether to save the plot as a PNG file. Default is `false`.
- `pdf::Bool=false`: Whether to save the plot as a PDF file. Default is `false`.
- `svg::Bool=false`: Whether to save the plot as an SVG file. Default is `false`.
- `filename::String="phaseDiagram"`: The filename for the saved plot. Default is "phaseDiagram".

# Returns
- `fig`: The matplotlib figure object.

# Example
```julia
julia>
```
"""
function plot2D(
    nums::T1,
    param_range1::T2,
    param_range2::T3;
    labels::Bool=true,
    disp::Bool=true,
    png::Bool=false,
    pdf::Bool=false,
    svg::Bool=false,
    filename::String="phaseDiagram",
) where {T1<:AbstractArray,T2<:AbstractVector,T3<:AbstractVector}
    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p_1")
        ax.set_ylabel(L"p_2")
    end

    if disp == true || png == true || pdf == true || svg == true
        im = ax.imshow(
            nums;
            cmap="jet",
            interpolation="none",
            origin="lower",
            extent=(param_range1[1], param_range1[end], param_range2[1], param_range2[end]),
            aspect="auto",
        )
        fig.colorbar(im; ax=ax, ticks=matplotlib.ticker.MaxNLocator(; integer=true))
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    return fig
end

@doc raw"""
    plot2D(result::NamedTuple; labels::Bool=true, disp::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String="phaseDiagram")
    
Plot a 2D phase diagram.

# Arguments
- `result::NamedTuple`: A named tuple containing the result data.
- `labels::Bool=true`: Whether to display axis labels.
- `disp::Bool=true`: Whether to display the plot.
- `png::Bool=false`: Whether to save the plot as a PNG file.
- `pdf::Bool=false`: Whether to save the plot as a PDF file.
- `svg::Bool=false`: Whether to save the plot as an SVG file.
- `filename::String="phaseDiagram"`: The filename for the saved plot.

# Returns
- `fig`: The matplotlib figure object.

# Example
```julia
julia>
```
"""
function plot2D(
    result::NamedTuple;
    labels::Bool=true,
    disp::Bool=true,
    png::Bool=false,
    pdf::Bool=false,
    svg::Bool=false,
    filename::String="phaseDiagram",
)
    fig = figure()
    ax = fig.add_subplot(111)

    if labels == true
        ax.set_xlabel(L"p_1")
        ax.set_ylabel(L"p_2")
    end

    nums_half = if result.nums isa Array{Float64,2}
        transpose(result.nums)
    else
        transpose(sum(@view(result.nums[1:(end ÷ 2), :, :]); dims=1)[1, :, :])
    end

    if disp == true || png == true || pdf == true || svg == true
        im = ax.imshow(
            nums_half;
            cmap="jet",
            interpolation="none",
            origin="lower",
            extent=(
                result.param1[1], result.param1[end], result.param2[1], result.param2[end]
            ),
            aspect="auto",
        )
        fig.colorbar(im; ax=ax, ticks=matplotlib.ticker.MaxNLocator(; integer=true))
    end

    p = (; disp, png, pdf, svg, filename)
    output(fig, p)
    return fig
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
        # plotshow()
        display(fig)
    end
end
