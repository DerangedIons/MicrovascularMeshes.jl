#--------------------------------------------------------------------------------# Scaling-law trees

"""
    ScalingTree{T}

A coordinate-free bifurcating vessel tree from [`scaling_tree`](@ref): topology, a generation
for each vessel, and its diameter and length. Nothing is placed in space; a caller embeds it
or hands it to a solver that only needs the hydraulics.

Vessels are numbered breadth-first, the inlet first: vessel `i` (`i > 1`) is fed by vessel
`i ÷ 2`, so generation `g` is the vessels `2^g:2^(g+1)-1`.

### Fields
- `parents::Vector{Int}`: the vessel that feeds each vessel, `0` for the inlet.
- `generation::Vector{Int}`: forks between the inlet and each vessel, `0` for the inlet.
- `diameter::Vector{T}`: the vessel's diameter, in `leaf_diameter`'s units.
- `length::Vector{T}`: the vessel's length, in `leaf_length`'s units.
"""
struct ScalingTree{T<:Real}
    parents::Vector{Int}
    generation::Vector{Int}
    diameter::Vector{T}
    length::Vector{T}
end

"""
    scaling_tree([rng=default_rng(),] depth; leaf_diameter, leaf_length, β, γ=3,
                 cv_diameter=0, cv_length=0,
                 diameter_distribution=nothing, length_distribution=nothing) -> ScalingTree

A symmetric tree of `depth` bifurcations whose vessel diameters and lengths follow measured
scaling laws: each generation has a mean diameter and a mean length set by two exponents, and
every vessel is drawn about those means with its own generation's spread.

**Why these statistics.** A vessel's Poiseuille resistance goes as `L/D⁴`, so the
per-generation statistics of `D` and `L` set where along the tree the pressure drops and how
unevenly flow divides among the leaves.

- **Diameter.** A fork conserves `D^γ`, `D_parent^γ = Σ D_daughter^γ`, so generation `g` (`0`
  the inlet) has mean `D̄_g = leaf_diameter · 2^((depth − g)/γ)`. Murray's `γ = 3` minimizes
  pumping plus blood-volume cost (Murray 1926), and arterioles below about 100 µm obey it:
  `Q ∝ D^3.01` (99.9 % CI 2.86–3.14) over 6–108 µm in rat cremaster (Mayrovitz & Roy 1983),
  and `x ≈ 3` in peripheral resistance arteries across studies (Nakamura & Awa 2014). The
  smaller values reported for large trees, `Q ∝ D^(7/3)` (Huo & Kassab 2012) or about `D²` in
  pig coronaries (Choy & Kassab 2008), regress each subtree's stem flow on its stem diameter
  across a whole tree, which is not a per-fork exponent. A lower `γ` widens each parent
  against its daughters and moves resistance toward the leaves.
- **Length.** `L̄_g = leaf_length · (D̄_g / leaf_diameter)^β`, and `β` has no default because
  its value depends on what one vessel stands for. A Strahler *element*, a run of segments of
  one order with its side branches, gives `β ≈ 1–1.2` (Van Bavel & Spaan 1992; Huo & Kassab
  2012), and about 2 in hamster skeletal-muscle arterioles under a dorsal skin-fold window
  (Lapi et al. 2019). A *segment* between two forks gives `β ≈ 0`: segment length barely
  tracks diameter (Wischgoll et al. 2009). A symmetric tree has no side branches, so read a
  generation as an element when the tree stands in for a deeper asymmetric one, and as a
  segment when it is the anatomy itself. With `γ = 3`, `β = 1` gives every generation the same
  pressure drop (in general `β = 4 − γ` does); `β = 0` concentrates the drop at the leaves.
- **Spread.** Each vessel's length and diameter are its generation's mean times a lognormal
  factor `exp(σ_g z − σ_g²/2)` with `σ_g² = log(1 + cv_g²)` and `z` standard normal from `rng`,
  which has mean one and coefficient of variation exactly `cv_g`. Lengths within an order are
  lognormal (Schröder et al. 1990; Kuhnle et al. 1993). Within one arteriolar network
  `cv_diameter ≈ 0.1` and `cv_length ≈ 0.2–0.4`; across networks, or by segment, lengths
  spread by 0.5–1 (Koller et al. 1987; Lapi et al. 2019; Wischgoll et al. 2009).
  `cv_diameter` and `cv_length` are each a number for every generation or a vector of
  `depth + 1`, root first. Lengths are drawn before diameters, vessel by vessel, and a quantity
  whose CV is zero everywhere draws nothing.
- **Other distributions.** To scatter by something other than a lognormal, pass
  `diameter_distribution` or `length_distribution` in place of the CV: each vessel's factor is
  then `rand(rng, d)` from its generation's `d`, any sampler such as a Distributions.jl
  distribution, or a vector of numbers to resample measured factors uniformly. Give one sampler
  for every generation or a vector of `depth + 1`, root first. The factor scales the
  generation's mean, so `d` should have mean one: a mean of `m` multiplies the whole generation
  by `m` and moves it off `γ` and `β`. Every draw must be finite and positive.

**Independent draws overstate heterogeneity.** Every vessel is drawn independently, so `γ`
holds generation by generation on average, not at every fork. Real trees size a vessel to the
subtree it feeds: a pig coronary tree reconstructed from per-order statistics gives a terminal
flow CV of 0.95 against 0.18–0.20 measured, and 0.22 once diameters follow the downstream tree
(Namani et al. 2018); converting measured lung diameter asymmetry through Poiseuille predicts
136 % flow dispersion against about 35 % measured (Glenny et al. 1991). Measured CVs fed in
here therefore bound the leaf-flow heterogeneity from above rather than calibrate it.

**Orders are not generations.** Morphometry is tabulated by diameter-defined Strahler order,
and a real tree has about three vessels of one order per vessel of the next, not two:
branching ratios of 2.76 (Jiang et al. 1994), 3.30 (Van Bavel & Spaan 1992) and 2.72 (Lapi et
al. 2019). A per-order diameter or length ratio is therefore not a per-generation one. Carry a
table over through its exponents (`γ`, `β`) and its per-order spreads matched by diameter,
never through its order-to-order ratios.

### Keyword Arguments
- `leaf_diameter`, `leaf_length`: the mean diameter and length of the last generation.
- `β`: the length–diameter exponent; required.
- `γ`: the junction exponent, positive.
- `cv_diameter`, `cv_length`: coefficient of variation per generation, `≥ 0`.
- `diameter_distribution`, `length_distribution`: the per-vessel factor's distribution, in
  place of the matching CV.

### Examples
```julia
scaling_tree(3; leaf_diameter=30.0, leaf_length=600.0, β=1)               # Murray, 20 diameters long
scaling_tree(6; leaf_diameter=30.0, leaf_length=600.0, β=1.2)             # element exponents
scaling_tree(Xoshiro(1), 8; leaf_diameter=30.0, leaf_length=600.0, β=1.2,
             cv_diameter=0.1, cv_length=0.3)                               # within-network spread
scaling_tree(Xoshiro(1), 8; leaf_diameter=30.0, leaf_length=600.0, β=1.2,
             length_distribution=Gamma(11, 1 / 11))                       # Gamma, CV 0.3 (Distributions.jl)
```

### References
- Choy JS, Kassab GS (2008) J Appl Physiol 104:1281. PMID 18323461
- Glenny RW, Robertson HT, Yamashiro S, Bassingthwaighte JB (1991) J Appl Physiol 70:2351. PMID 1885430
- Huo Y, Kassab GS (2012) J R Soc Interface 9:190. PMID 21676970
- Jiang ZL, Kassab GS, Fung YC (1994) J Appl Physiol 76:882. PMID 8175603
- Koller A, Dawant B, Liu A, Popel AS, Johnson PC (1987) Am J Physiol 253:H154. PMID 3605362
- Kuhnle GE et al. (1993) Int J Microcirc Clin Exp 12:313. PMID 8375965
- Lapi D et al. (2019) Front Physiol 9:1953. PMID 30713505
- Mayrovitz HN, Roy J (1983) Am J Physiol 245:H1031. PMID 6660303
- Murray CD (1926) Proc Natl Acad Sci USA 12:207. PMID 16576980
- Nakamura Y, Awa S (2014) Physiol Rep 2:e00236. PMID 24744905
- Namani R, Kassab GS, Lanir Y (2018) Front Physiol 9:1069. PMID 30210353
- Schröder S et al. (1990) Fortschr Ophthalmol 87:52. PMID 2323698
- Van Bavel E, Spaan JA (1992) Circ Res 71:1200. PMID 1394880
- Wischgoll T, Choy JS, Kassab GS (2009) Am J Physiol Heart Circ Physiol 297:H1949. PMID 19749169
"""
function scaling_tree(
    rng::AbstractRNG,
    depth::Integer;
    leaf_diameter::Real,
    leaf_length::Real,
    β::Real,
    γ::Real=3,
    cv_diameter=0,
    cv_length=0,
    diameter_distribution=nothing,
    length_distribution=nothing,
)
    depth ≥ 0 || throw(ArgumentError("depth must be ≥ 0, got $depth"))
    isfinite(γ) && γ > 0 || throw(ArgumentError("γ = $γ is the junction exponent and must be positive and finite"))
    isfinite(β) || throw(ArgumentError("β = $β is the length–diameter exponent and must be finite"))
    T = float(promote_type(typeof(leaf_diameter), typeof(leaf_length)))
    n = 2^(depth + 1) - 1
    parents = [0; (2:n) .÷ 2]
    generation = map(i -> ndigits(i; base=2) - 1, 1:n)
    (D_leaf, L_leaf) = (T(leaf_diameter), T(leaf_length))
    mean_diameter = map(g -> D_leaf * T(2)^((depth - g) / T(γ)), generation)
    mean_length = map(d -> L_leaf * (d / D_leaf)^T(β), mean_diameter)
    len = mean_length .* _factors(rng, T, cv_length, length_distribution, depth, generation, :length)
    diameter = mean_diameter .* _factors(rng, T, cv_diameter, diameter_distribution, depth, generation, :diameter)
    return ScalingTree{T}(parents, generation, diameter, len)
end
scaling_tree(depth::Integer; kwargs...) = scaling_tree(default_rng(), depth; kwargs...)

# One factor per vessel about its generation's mean: a lognormal of the given CV, or a draw from
# the given distribution when there is one.
function _factors(rng::AbstractRNG, ::Type{T}, cv, ::Nothing, depth, generation, name) where {T}
    return _lognormal_factors(rng, _per_generation(T, cv, depth, Symbol(:cv_, name)), generation)
end
function _factors(rng::AbstractRNG, ::Type{T}, cv, distribution, depth, generation, name) where {T}
    all(iszero, cv) || throw(ArgumentError("give cv_$name or $(name)_distribution, not both"))
    d = _per_generation(distribution, depth, Symbol(name, :_distribution))
    return map(generation) do g
        x = T(rand(rng, d[g+1]))
        isfinite(x) && x > 0 || throw(ArgumentError("$(name)_distribution drew $x; a factor must be finite and positive"))
        x
    end
end

# One mean-one lognormal factor per vessel, whose coefficient of variation is exactly its
# generation's `cv[g+1]`. A CV of zero everywhere draws nothing from `rng`.
function _lognormal_factors(rng::AbstractRNG, cv::Vector{T}, generation::Vector{Int}) where {T}
    all(c -> isfinite(c) && c ≥ 0, cv) || throw(ArgumentError("a coefficient of variation must be finite and ≥ 0, got $cv"))
    all(iszero, cv) && return ones(T, length(generation))
    σ = map(c -> sqrt(log1p(c^2)), cv)
    return map(g -> exp(σ[g+1] * randn(rng, T) - σ[g+1]^2 / 2), generation)
end

# A per-generation quantity as a vector of `depth + 1`, root first.
_per_generation(::Type{T}, x::Real, depth, name) where {T} = fill(T(x), depth + 1)
_per_generation(::Type{T}, x::AbstractVector{<:Real}, depth, name) where {T} = collect(T, _check_generations(x, depth, name))

# A per-generation distribution likewise. A vector of numbers is one empirical distribution, not
# one number per generation.
_per_generation(d, depth, name) = fill(d, depth + 1)
_per_generation(d::AbstractVector, depth, name) = eltype(d) <: Real ? fill(d, depth + 1) : _check_generations(d, depth, name)

function _check_generations(x::AbstractVector, depth, name)
    length(x) == depth + 1 || throw(
        DimensionMismatch("a tree of depth $depth has $(depth + 1) generations, but $name has $(length(x)) entries"),
    )
    return x
end
