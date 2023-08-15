var documenterSearchIndex = {"docs":
[{"location":"2D/flux/#Two-dimensional-square-lattice-model-with-flux","page":"Square Lattice w/ Flux","title":"Two-dimensional square lattice model with flux","text":"","category":"section"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"A two-dimensional example is presented here:","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"julia> function H₀(k, p) # landau\n    k1, k2 = k\n    t = 1\n\n    Hsize = 6\n    Hmat = zeros(ComplexF64, Hsize, Hsize)\n\n    ϕ = 2π * p / Hsize\n\n    for i in 1:Hsize\n        Hmat[i, i] = -2t * cos(k2 - i * ϕ)\n    end\n\n    for i in 1:Hsize-1\n        Hmat[i, i+1] = -t\n        Hmat[i+1, i] = -t\n    end\n\n    Hmat[1, Hsize] = -t * exp(-im * Hsize * k1)\n    Hmat[Hsize, 1] = -t * exp(im * Hsize * k1)\n\n    Hmat\nend","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"To calculate the dispersion, run:","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"julia> H(k) = H₀(k, 1)\njulia> showBand(H; value=false, disp=true)","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"(Image: Dispersion of 2D square lattice with flux model)","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"Then we can compute the Chern numbers using calcChern:","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"julia> calcChern(H)","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"The output is:","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"(TopologicalNumber = [6, 6, -12, -12, 6, 6], Total = 0)","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the first Chern number for each band.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the first Chern numbers for each band. Total is a quantity that should always return zero.","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"One-dimensional phase diagram is given by:","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"julia> param = 1:6\njulia> calcPhaseDiagram(H₀, param, \"Chern\"; plot=true)","category":"page"},{"location":"2D/flux/","page":"Square Lattice w/ Flux","title":"Square Lattice w/ Flux","text":"%%%dummy%%% (Image: One-dimensional phase diagram) %%%dummy%%%","category":"page"},{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"","category":"page"},{"location":"1D/SSH/#The-Su-Schriffer-Heeger-(SSH)-model","page":"SSH model","title":"The Su-Schriffer-Heeger (SSH) model","text":"","category":"section"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"Here's a simple example of the SSH Hamiltonian:","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"julia> using TopologicalNumbers\njulia> function H₀(k, p)\n            [\n                0 p[1]+p[2]*exp(-im * k)\n                p[1]+p[2]*exp(im * k) 0\n            ]\n        end","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"The band structure is computed as follows:","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"julia> H(k) = H₀(k, (0.9, 1.0))\njulia> showBand(H; value=false, disp=true)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"(Image: Band structure of SSH model)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"Next, we can calculate the winding numbers using calcBerryPhase:","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"julia> calcBerryPhase(H)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"The output is:","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"(TopologicalNumber = [1, 1], Total = 0)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the winding number for each band.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the winding numbers for each band (mod 2). Total is a quantity that should always return zero.","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"One-dimensional phase diagram is given by:","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"julia> H(k, p) = H₀(k, (p, 1.0))\n\njulia> param = range(-2.0, 2.0, length=1001)\njulia> calcPhaseDiagram(H, param, \"BerryPhase\"; plot=true)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"(Image: One-dimensional phase diagram of SSH model)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"Also, two-dimensional phase diagram is given by:","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"julia> param = range(-2.0, 2.0, length=101)\njulia> calcPhaseDiagram(H₀, param, param, \"BerryPhase\"; plot=true)","category":"page"},{"location":"1D/SSH/","page":"SSH model","title":"SSH model","text":"(Image: Two-dimensional phase diagram of SSH model)","category":"page"},{"location":"lib/public/#Public","page":"Public","title":"Public","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"","category":"page"},{"location":"lib/public/","page":"Public","title":"Public","text":"Modules = [TopologicalNumbers]\nInternal = false","category":"page"},{"location":"lib/public/#TopologicalNumbers.calcBerryPhase-Tuple{Function}","page":"Public","title":"TopologicalNumbers.calcBerryPhase","text":"Calculate the winding numbers in the one-dimensional case.\n\ncalcBerryPhase(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)\n\nArguments\n\nHamiltonian::Function: the Hamiltonian matrix function with one-dimensional wavenumber k as an argument.\nN::Int: the number of meshes when discretizing the Brillouin Zone. It is preferable for N to be an odd number to increase the accuracy of the calculation.\ngapless::Real: the threshold that determines the state to be degenerate. Coarsening the mesh(N) but increasing gapless will increase the accuracy of the calculation.\nrounds::Bool: an option to round the value of the topological number to an integer value. The topological number returns a value of type Int when true, and a value of type Float when false.\n\nDefinition\n\nThe Berry phase of the nth band nu_n is defined by\n\nnu_n=frac1pi iint_mathrmBZdkA_n(k)\n\nThe integral range mathrmBZ(Brillouin Zone) is kin02pi. A_n(k) is the Berry conection at wavenumber k.\n\nA_n(k)=braPsi_n(k)partial_kketPsi_n(k)\n\nketPsi_n(k) is the wave function of the nth band.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.calcChern-Tuple{Function}","page":"Public","title":"TopologicalNumbers.calcChern","text":"Calculate the first Chern numbers in the two-dimensional case with reference to Fukui-Hatsugai-Suzuki method [1].\n\ncalcChern(Hamiltonian::Function; N::Int=51, gapless::Real=0.0, rounds::Bool=true)\n\nArguments\n\nHamiltionian::Function: the Hamiltonian matrix with one-dimensional wavenumber k as an argument.\nN::Int=51: The number of meshes when discretizing the Brillouin Zone. It is preferable for N to be an odd number to increase the accuracy of the calculation.\ngapless::Real: The threshold that determines the state to be degenerate. Coarsening the mesh(N) but increasing gapless will increase the accuracy of the calculation.\nrounds::Bool=true: An option to round the value of the topological number to an integer value. The topological number returns a value of type Int when true, and a value of type Float when false.\n\nDefinition\n\nThe firs Chern number of the nth band nu_n is defined by\n\nnu_n=frac12pi iint_mathrmBZdbmkleft(partial_k_1A_n2(bmk)-partial_k_2A_n1(bmk)right)\n\nThe integral range mathrmBZ(Brillouin Zone) is bmkin02pi^2. A_ni(bmk) is the Berry connection at wavenumber bmk.\n\nA_ni(bmk)=braPsi_n(bmk)partial_k_iketPsi_n(bmk)\n\nketPsi_n(bmk) is the wave function of the nth band.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.calcPhaseDiagram-Union{Tuple{T}, Tuple{Function, T, String}} where T<:(AbstractVector)","page":"Public","title":"TopologicalNumbers.calcPhaseDiagram","text":"calcPhaseDiagram(H::Function, param_range::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.calcPhaseDiagram-Union{Tuple{T}, Tuple{Function, T, T, String}} where T<:(AbstractVector)","page":"Public","title":"TopologicalNumbers.calcPhaseDiagram","text":"calcPhaseDiagram(H::Function, param_range1::T, param_range2::T, alg::String; N::Int=51, gapless::Real=0.0, rounds::Bool=true, plot::Bool=false) where {T<:AbstractVector}\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.calcZ2-Tuple{Function}","page":"Public","title":"TopologicalNumbers.calcZ2","text":"Calculate the mathbbZ_2 numbers in the two-dimensional case with reference to Shiozaki method [2].\n\ncalcZ2(Hamiltonian::Function; N::Int=50, rounds::Bool=true, TR::Bool=false)\n\nArguments\n\nHamiltonian::Function is a matrix with one-dimensional wavenumber k as an argument.\nN::Int is the number of meshes when discretizing the Brillouin Zone. It is preferable for N to be an odd number to increase the accuracy of the calculation.\nrounds::Bool is an option to round the value of the topological number to an integer value. The topological number returns a value of type Int when true, and a value of type Float when false.\n\nDefinition\n\nThe mathbbZ_2 number of the 2nth (and 2n-1th) band nu_n is defined by\n\nnu_n=frac12pi iint_mathrmBZdbmkleft(partial_k_1A_n2(bmk)-partial_k_2A_n1(bmk)right)-frac12pi ioint_partialmathrmBZdk_iA_njneq i(bmk)\n\nThe integral range mathrmBZ is bmkin02pitimes0pi half of BZ(Brillouin Zone). A_ni(bmk) is the Berry connection at wavenumber bmk.\n\nA_ni(bmk)=braPsi_n(bmk)partial_k_iketPsi_n(bmk)\n\nketPsi_n(bmk) is the wave function of the 2nth (and 2n-1th) band.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.plot1D-Tuple{NamedTuple}","page":"Public","title":"TopologicalNumbers.plot1D","text":"plot1D(result::NamedTuple; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String=\"phaseDiagram\")\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.plot1D-Union{Tuple{T}, Tuple{Matrix, T}} where T<:(AbstractVector)","page":"Public","title":"TopologicalNumbers.plot1D","text":"plot1D(nums::Matrix, param_range::T; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String=\"phaseDiagram\") where {T<:AbstractVector}\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.plot2D-Tuple{NamedTuple}","page":"Public","title":"TopologicalNumbers.plot2D","text":"plot2D(result::NamedTuple; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String=\"phaseDiagram\")\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.plot2D-Union{Tuple{T2}, Tuple{T1}, Tuple{T1, T2, T2}} where {T1<:AbstractArray, T2<:(AbstractVector)}","page":"Public","title":"TopologicalNumbers.plot2D","text":"plot2D(nums::T1, param_range1::T2, param_range2::T2; labels::Bool=true, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String=\"phaseDiagram\") where {T1<:AbstractArray,T2<:AbstractVector}\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#TopologicalNumbers.showBand-Tuple{Function}","page":"Public","title":"TopologicalNumbers.showBand","text":"Drawing the band structure of the Hamiltonian.\n\nshowBand(Hamiltonian::Function; N::Int=51, labels::Bool=true, value::Bool=true, disp::Bool=false, png::Bool=false, pdf::Bool=false, svg::Bool=false, filename::String=\"Band\")\n\nArguments\n\nHamiltonian::Function: the Hamiltonian matrix function of wave number bm k.\nN::Int: the number of divisions in the wave number space.\nlabels::Bool: whether to display the labels of the figure.\nvalue::Bool: whether to output the values of the wave number and the energy in the matrix form.\ndisp::Bool: whether to display the figure.\npng::Bool: whether to save the figure as a PNG file.\npdf::Bool: whether to save the figure as a PDF file.\nsvg::Bool: whether to save the figure as a SVG file.\nfilename::String: the name of the output file.\n\n\n\n\n\n\n\n","category":"method"},{"location":"examples_old/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples_old/#One-dimensional-case","page":"Examples","title":"One-dimensional case","text":"","category":"section"},{"location":"examples_old/#The-Su-Schriffer-Heeger-(SSH)-model","page":"Examples","title":"The Su-Schriffer-Heeger (SSH) model","text":"","category":"section"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"Here's a simple example of the SSH Hamiltonian:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> using TopologicalNumbers\njulia> function H(k) # set SSH Hamiltonian function of wavenumber k\n    g = 0.9\n    \n    [\n        0 g+exp(-im*k)\n        g+exp(im*k) 0\n    ]\nend","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The band structure is computed as follows:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> showBand(H; value=false, disp=true)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(Image: Band structure of SSH model)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"In this case, 1 signifies the dimension of the wavenumber space.","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"Next, we can calculate the winding numbers using calcBerryPhase:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> calcBerryPhase(H)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The output is:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(TopologicalNumber = [1, 1], Total = 0)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the winding number for each band.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the winding numbers for each band (mod 2). Total is a quantity that should always return zero.","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"One-dimensional phase diagram is given by:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> function H0(k, p)\n            [\n                0 p[1]+p[2]*exp(-im * k)\n                p[1]+p[2]*exp(im * k) 0\n            ]\n        end\njulia> H(k, p) = H0(k, (p, 1.0))\n\njulia> param = range(-2.0, 2.0, length=1001)\njulia> calcPhaseDiagram(H, param, \"BerryPhase\"; plot=true)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(Image: One-dimensional phase diagram of SSH model)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"Also, two-dimensional phase diagram is given by:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> param = range(-2.0, 2.0, length=101)\njulia> calcPhaseDiagram(H0, param, param, \"BerryPhase\"; plot=true)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(Image: Two-dimensional phase diagram of SSH model)","category":"page"},{"location":"examples_old/#Chern-numbers","page":"Examples","title":"Chern numbers","text":"","category":"section"},{"location":"examples_old/#Two-dimensional-square-lattice-with-flux-model","page":"Examples","title":"Two-dimensional square lattice with flux model","text":"","category":"section"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"A two-dimensional example is presented here:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> function H(k) # landau\n    k1, k2 = k\n    t = 1\n\n    Hsize = 6\n    Hmat = zeros(ComplexF64, Hsize, Hsize)\n\n    for i in 1:Hsize\n        Hmat[i, i] = -2*cos(k2-2pi*i/Hsize)\n    end\n\n    for i in 1:Hsize-1\n        Hmat[i, i+1] = -t\n        Hmat[i+1, i] = -t\n    end\n\n    Hmat[1, Hsize] = -t*exp(-im*k1)\n    Hmat[Hsize, 1] = -t*exp(im*k1)\n    \n    Hmat\nend","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"To calculate the dispersion, run:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> showBand(H; value=false, disp=true)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(Image: Dispersion of 2D square lattice with flux model)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"Then we can compute the Chern numbers using calcChern:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> calcChern(H)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The output is:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(TopologicalNumber = [1, 1, -2, -2, 1, 1], Total = 0)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the first Chern number for each band.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the first Chern numbers for each band. Total is a quantity that should always return zero.","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"One-dimensional phase diagram is given by:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> function H(k, p)\n    k1, k2 = k\n    t = p\n\n    Hsize = 6\n    Hmat = zeros(ComplexF64, Hsize, Hsize)\n\n    for i in 1:Hsize\n        Hmat[i, i] = -2 * cos(k2 - 2pi * i / Hsize)\n    end\n\n    for i in 1:Hsize-1\n        Hmat[i, i+1] = -t\n        Hmat[i+1, i] = -t\n    end\n\n    Hmat[1, Hsize] = -t * exp(-im * k1)\n    Hmat[Hsize, 1] = -t * exp(im * k1)\n\n    Hmat\nend\n\njulia> param = range(-2.0, 2.0, length=500)\njulia> calcPhaseDiagram(H, param, \"Chern\"; plot=true)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(Image: One-dimensional phase diagram)","category":"page"},{"location":"examples_old/#\\mathbb{Z}_2-numbers","page":"Examples","title":"mathbbZ_2 numbers","text":"","category":"section"},{"location":"examples_old/#The-Bernevig-Hughes-Zhang-(BHZ)-model","page":"Examples","title":"The Bernevig-Hughes-Zhang (BHZ) model","text":"","category":"section"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"As an example of a two-dimensional topological insulator, the BHZ model is presented here:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> function H(k) # BHZ\n    k1, k2 = k\n    t = 1\n\n    R0 = -2(cos(k1) + cos(k2)) + 1\n    R3 = 2sin(k2)\n    R4 = 2sin(k1)\n    R5 = -2t*(cos(k1) + cos(k2)) + 1\n\n    s0 = [1 0; 0 1]\n    sx = [0 1; 1 0]\n    sy = [0 -im; im 0]\n    sz = [1 0; 0 -1]\n\n    a0 = kron(s0, s0)\n    a1 = kron(sx, sx)\n    a2 = kron(sx, sy)\n    a3 = kron(sx, sz)\n    a4 = kron(sy, s0)\n    a5 = kron(sz, s0)\n\n    R0*a0+R3*a3+R4*a4+R5*a5\nend","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"To calculate the dispersion, execute:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> showBand(H; value=false, disp=true)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(Image: Dispersion of BHZ model)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"Next, we can compute the mathbbZ_2 numbers using calcZ2:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"julia> calcZ2(H)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The output is:","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"(TopologicalNumber = [1, 1], Total = 0)","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the mathbbZ_2 number for each each pair of two energy bands.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the mathbbZ_2 numbers for each pair of two energy bands. Total is a quantity that should always return zero.","category":"page"},{"location":"examples_old/","page":"Examples","title":"Examples","text":"Total is a value that should consistently return zero.","category":"page"},{"location":"2D/Haldane/#Haldane-model","page":"Haldane model","title":"Haldane model","text":"","category":"section"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"Hamiltonian of Haldane model is given by:","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"julia> function H₀(k, p) # landau\n           k1, k2 = k\n           J = 1.0\n           K = 1.0\n           ϕ, M = p\n\n           h0 = 2K * cos(ϕ) * (cos(k1) + cos(k2) + cos(k1 + k2))\n           hx = J * (1 + cos(k1) + cos(k2))\n           hy = J * (-sin(k1) + sin(k2))\n           hz = M - 2K * sin(ϕ) * (sin(k1) + sin(k2) - sin(k1 + k2))\n\n           s0 = [1 0; 0 1]\n           sx = [0 1; 1 0]\n           sy = [0 -im; im 0]\n           sz = [1 0; 0 -1]\n\n           h0 .* s0 .+ hx .* sx .+ hy .* sy .+ hz .* sz\n       end","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"The band structure is computed as follows:","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"julia> H(k) = H₀(k, (π/3, 0.5))\njulia> showBand(H; value=false, disp=true)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"(Image: Band structure of Haldane model)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"Then we can compute the Chern numbers using calcChern:","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"julia> calcChern(H)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"The output is:","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"(TopologicalNumber = [1, -1], Total = 0)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the first Chern number for each band.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the first Chern numbers for each band. Total is a quantity that should always return zero.","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"One-dimensional phase diagram is given by:","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"julia> H(k, p) = H₀(k, (p, 2.5))\n\njulia> param = range(-π, π, length=1000)\njulia> calcPhaseDiagram(H, param, \"Chern\"; plot=true)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"(Image: One-dimensional phase diagram of Haldane model)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"Also, two-dimensional phase diagram is given by:","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"julia> param1 = range(-π, π, length=100)\njulia> param2 = range(-6.0, 6.0, length=100)\njulia> calcPhaseDiagram(H₀, param1, param2, \"Chern\"; plot=true)","category":"page"},{"location":"2D/Haldane/","page":"Haldane model","title":"Haldane model","text":"(Image: Two-dimensional phase diagram of Haldane model)","category":"page"},{"location":"lib/internal/#Internal","page":"Internal","title":"Internal","text":"","category":"section"},{"location":"lib/internal/","page":"Internal","title":"Internal","text":"Modules = [TopologicalNumbers]\nPublic = false","category":"page"},{"location":"2D/BHZ/#The-Bernevig-Hughes-Zhang-(BHZ)-model","page":"BHZ model","title":"The Bernevig-Hughes-Zhang (BHZ) model","text":"","category":"section"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"As an example of a two-dimensional topological insulator, the BHZ model is presented here:","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"julia> function H₀(k, p) # BHZ\n    k1, k2 = k\n    tₛₚ = 1\n    t₁ = ϵ₁ = 2\n    ϵ₂, t₂ = p\n\n    R0 = -t₁*(cos(k1) + cos(k2)) + ϵ₁/2\n    R3 = 2tₛₚ*sin(k2)\n    R4 = 2tₛₚ*sin(k1)\n    R5 = -t₂*(cos(k1) + cos(k2)) + ϵ₂/2\n\n    s0 = [1 0; 0 1]\n    sx = [0 1; 1 0]\n    sy = [0 -im; im 0]\n    sz = [1 0; 0 -1]\n\n    a0 = kron(s0, s0)\n    a1 = kron(sx, sx)\n    a2 = kron(sx, sy)\n    a3 = kron(sx, sz)\n    a4 = kron(sy, s0)\n    a5 = kron(sz, s0)\n\n    R0 .* a0 .+ R3 .* a3 .+ R4 .* a4 .+ R5 .* a5\nend","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"To calculate the dispersion, execute:","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"julia> H(k) = H₀(k, (2, 2))\njulia> showBand(H; value=false, disp=true)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"(Image: Dispersion of BHZ model)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"Next, we can compute the mathbbZ_2 numbers using calcZ2:","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"julia> calcZ2(H)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"The output is:","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"(TopologicalNumber = [1, 1], Total = 0)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"The first argument TopologicalNumber in the named tuple is an vector that stores the mathbbZ_2 number for each pair of two energy bands.  The vector is arranged in order of bands, starting from the one with the lowest energy. The second argument Total stores the total of the mathbbZ_2 numbers for each pair of two energy bands. Total is a quantity that should always return zero.","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"One-dimensional phase diagram is given by:","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"julia> H(k, p) = H₀(k, (p, 0.25))\n\njulia> param = range(-2, 2, length=1000)\njulia> calcPhaseDiagram(H, param, \"Z2\"; plot=true)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"(Image: One-dimensional phase diagram of BHZ model)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"Also, two-dimensional phase diagram is given by:","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"julia> param1 = range(-2, 2, length=100)\njulia> param2 = range(-0.5, 0.5, length=100)\njulia> calcPhaseDiagram(H₀, param1, param2, \"Z2\"; plot=true)","category":"page"},{"location":"2D/BHZ/","page":"BHZ model","title":"BHZ model","text":"(Image: Two-dimensional phase diagram of BHZ model)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = TopologicalNumbers","category":"page"},{"location":"#TopologicalNumbers.jl","page":"Home","title":"TopologicalNumbers.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TopologicalNumbers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage)","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TopologicalNumbers.jl is a Julia package designed to calculate topological numbers, such as the Chern numbers and mathbbZ_2 numbers,  using a numerical approach based on the Fukui-Hatsugai-Suzuki method [1] or the Shiozaki method [2].","category":"page"},{"location":"","page":"Home","title":"Home","text":"This software is released under the MIT License, please see the LICENSE file for more details.   It is confirmed to work on Julia 1.6 (LTS) and 1.9.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install TopologicalNumbers.jl, run the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add TopologicalNumbers","category":"page"},{"location":"","page":"Home","title":"Home","text":"Alternatively, you can use:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"TopologicalNumbers\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nIf you are using a headless server, you may receive an error related to the GLMakie package. To resolve this issue, please refer to the Makie documentation or the GLMakie troubleshooting guide.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package includes the following functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"showBand to calculate the dispersion relation,\ncalcBerryPhase to calculate the winding numbers in the one-dimensional case,\ncalcChern to calculate the first Chern numbers in the two-dimensional case,\ncalcZ2 to calculate the mathbbZ_2 numbers in the two-dimensional case,\ncalcPhaseDiagram to calculate the phase diagram using the several methods.","category":"page"}]
}
