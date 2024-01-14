@doc raw"""
"""
abstract type TopologicalNumbersAlgorithms end

# Algorithms for calculating the Berry phase
@doc raw"""
"""
abstract type BerryPhaseAlgorithms <: TopologicalNumbersAlgorithms end
# struct Int1DBP <: Z2Algorithms end

@doc raw"""
Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite)

# Definition

The Berry phase of the $n$th band $\nu_{n}$ is defined by
```math
\nu_{n}=\frac{1}{\pi}\sum_{k\in\mathrm{BZ}}U_{n}(k)
```
The range $\mathrm{BZ}$(Brillouin Zone) is $k\in[0,2\pi]$. $U_{n,i}(k)$ is the link variable at wavenumber $k$. $e_{1}$ is the unit vector.
```math
U_{n}(k)=\braket{\Psi_{n}(k)|\Psi_{n}(k+e_{1})}
```
$\ket{\Psi_{n}(k)}$ is the wave function of the $n$th band.
"""
struct BP <: BerryPhaseAlgorithms end

# Algorithms for calculating the first Chern number
@doc raw"""
"""
abstract type FirstChernAlgorithms <: TopologicalNumbersAlgorithms end
# struct IntFChern <: FirstChernAlgorithms end

@doc raw"""
Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite)

# Definition
 The first Chern number of the $n$th band $\nu_{n}$ is defined by
```math
\nu_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,2}^{*}(\bm{k})\right]\right]
```
 The range $\mathrm{BZ}$(Brillouin Zone) is $\bm{k}\in[0,2\pi]^{2}$. $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
"""
struct FHS <: FirstChernAlgorithms end


# Algorithms for calculating the second Chern number
@doc raw"""
"""
abstract type SecondChernAlgorithms <: TopologicalNumbersAlgorithms end
# struct IntSChern <: SecondChernAlgorithms end

@doc raw"""
Fukui-Hatsugai-Suzuki method [Mochol-Grzelak2018Efficient](@cite)
"""
struct FHS2 <: SecondChernAlgorithms end

# Algorithms for calculating the Z2 invariant
@doc raw"""
"""
abstract type Z2Algorithms <: TopologicalNumbersAlgorithms end
# struct Int2DZ2 <: Z2Algorithms end

@doc raw"""
Shiozaki method [Shiozaki2023discrete](@cite)

# Definition
 The $\mathbb{Z}_{2}$ number of the $2n$th (and $2n-1$th) band $\nu_{n}$ is defined by
```math
\nu_{n}=F_{n}-\left(P_{n}(0)-P_{n}(\pi)\right)
```
 $F_{n}$ is the Berry flux of the $n$th band in the $\mathrm{BZ}'$. The range $\mathrm{BZ}'$ is $\bm{k}\in[0,2\pi]\times[0,\pi]$ half of BZ(Brillouin Zone).
```math
F_{n}=\frac{1}{2\pi}\sum_{\bm{k}\in\mathrm{BZ}'}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
 $P_{n}(k_{2})$ is the time-reversal polarization at wavenumber $k_{2}$.
```math
P_{n}(k_{2})=\frac{1}{2\pi}\frac{\mathrm{Pf}[\omega(0,k_{2})]}{\mathrm{Pf}[\omega(\pi,k_{2})]}\sum_{k_{1}=0}^{\pi-e_{1}}U_{n,1}(\bm{k})
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $2n$th (and $2n-1$th) band. $\omega(\bm{k})$ is the unitary matrix given by
```math
\omega(\bm{k})=\bra{\Psi(-\bm{k})}T\ket{\Psi(\bm{k})}
```
 $T$ is the time-reversal operator.
"""
struct Shio <: Z2Algorithms end


# Algorithms for calculating the local Berry flux
@doc raw"""
"""
abstract type BerryFluxAlgorithms <: TopologicalNumbersAlgorithms end

@doc raw"""
Fukui-Hatsugai-Suzuki method [Fukui2005Chern](@cite)

# Definition
 The Berry flux at the wavenumber $\bm{k}$ of the $n$th band $F_{n}(\bm{k})$ is defined by
```math
F_{n}(\bm{k})=f_{n}(\bm{k})-df_{n}(\bm{k})
```
```math
f_{n}(\bm{k})=\frac{1}{2\pi}\mathrm{Im}\left[\mathrm{Log}\left[U_{n,1}(\bm{k})U_{n,2}(\bm{k}+\bm{e}_{1})U_{n,1}^{*}(\bm{k}+\bm{e}_{2})U_{n,1}^{*}(\bm{k})\right]\right]
```
```math
df_{n}(\bm{k})=\frac{1}{2\pi}\mathrm{Im}\left[\mathrm{Log}[U_{n,1}(\bm{k})]+\mathrm{Log}[U_{n,2}(\bm{k}+\bm{e}_{1})]-\mathrm{Log}[U_{n,1}(\bm{k}+\bm{e}_{2})]-\mathrm{Log}[U_{n,1}(\bm{k})]\right]
```
 $U_{n,i}(\bm{k})$ is the link variable at wavenumber $\bm{k}$. $\bm{e}_{i}$ is the unit vector.
```math
U_{n,i}(\bm{k})=\braket{\Psi_{n}(\bm{k})|\Psi_{n}(\bm{k}+\bm{e}_{i})}
```
 $\ket{\Psi_{n}(\bm{k})}$ is the wave function of the $n$th band.
"""
struct FHSlocal2 <: BerryFluxAlgorithms end


# Algorithms for finding and calculating the Weyl points
@doc raw"""
"""
abstract type WeylPointsAlgorithms <: TopologicalNumbersAlgorithms end

@doc raw"""
"""
struct FHSsurface <: WeylPointsAlgorithms end

@doc raw"""
"""
struct FHSlocal3 <: WeylPointsAlgorithms end

@doc raw"""
"""
struct Evar <: WeylPointsAlgorithms end