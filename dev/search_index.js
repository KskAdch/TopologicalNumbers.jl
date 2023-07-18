var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TopologicalNumbers","category":"page"},{"location":"#TopologicalNumbers.jl","page":"Home","title":"TopologicalNumbers.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TopologicalNumbers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TopologicalNumbers]","category":"page"},{"location":"#TopologicalNumbers.jl:-A-Julia-package-for-calculating-topological-numbers","page":"Home","title":"TopologicalNumbers.jl: A Julia package for calculating topological numbers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage)","category":"page"},{"location":"#Abstract","page":"Home","title":"Abstract","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TopologicalNumbers.jl is a Julia package to calculate topological numbers, including the Chern numbers and mathbbZ_2 numbers using a numerical approach based on the Fukui-Hataugai-Suzuki method. This package serves the following functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Dispersion to calculate the dispersion relation,\nQuantizedBerryPhase to calculate the winding numbers in the one-dimensional case,\nFirstChern to calculate the first Chern numbers in the two-dimensional case,\nZ2invariants2D to calculate the mathbbZ_2 numbers in the two-dimensional case.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This software is released under the MIT License, see LICENSE. We checked that it works on Julia 1.6 (LTS) and 1.9.","category":"page"},{"location":"#Install","page":"Home","title":"Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install TopologicalNumbers.jl with:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"https://github.com/KskAdch/TopologicalNumbers.jl\")","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"#The-Su-Schriffer-Heeger-(SSH)-model","page":"Home","title":"The Su-Schriffer-Heeger (SSH) model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Here is a simple example of the SSH Hamiltonian:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using TopologicalNumbers\njulia> function H(k) # set SSH Hamiltonian function of wavenumber k\n    g = 0.9\n    \n    [\n        0 g+exp(-im*k)\n        g+exp(im*k) 0\n    ]\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"A band structure is here:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> Dispersion(H, 1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Band structure of SSH model)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here, 1 is the dimension of the wavenumber space.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then we can calculate the winding numbers using QuantizedBerryPhase:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> QuantizedBerryPhase(H)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Output is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(TopologicalNumber = [1, 1], Total = 0)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This means...(edit required)","category":"page"},{"location":"#Two-dimensional-square-lattice-with-flux-model","page":"Home","title":"Two-dimensional square lattice with flux model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Two-dimensional example is here:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> function H(k) # landau\n    k1, k2 = k\n    t = 1\n\n    Hsize = 6\n    Hmat = zeros(ComplexF64, Hsize, Hsize)\n\n    for i in 1:Hsize\n        Hmat[i, i] = -2*cos(k2-2pi*i/Hsize)\n    end\n\n    for i in 1:Hsize-1\n        Hmat[i, i+1] = -t\n        Hmat[i+1, i] = -t\n    end\n\n    Hmat[1, Hsize] = -t*exp(-im*k1)\n    Hmat[Hsize, 1] = -t*exp(im*k1)\n    \n    Hmat\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"To calcurate Dispersion, run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> Dispersion(H, 2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Dispersion of 2D square lattice with flux model)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then we can calculate the Chern numbers using FirstChern:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> FirstChern(H)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Output is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This means...(edit required)","category":"page"}]
}