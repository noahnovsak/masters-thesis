module ppt2

using Random
using LinearAlgebra
using JuMP
using MosekTools

using DynamicPolynomials: AbstractPolynomial, @polyvar, variables, coefficient
using SumOfSquares: SOSModel, SOSCone, gram_matrix

import Ket: partial_transpose, partial_trace, entanglement_robustness   # `using Ket` would clash with symmetric_projector

export pncp_mat, ampliation, rand_ppt, rand_sep, rand_psd, is_ppt,
    antisymmetric_projector, gram_freedom, is_block_positive,
    detect_trace, detect_ampliation, detect_dps, test_ppt2, min_ppt_witness,
    min_eig, has_negative_eig, has_negative_eig!, min_ppt2_witness,
    load_batches, load_meta, load_states, batch_id_of,
    completed_batches, batch_counts, write_meta!, sample_batch, generate_dataset


⊗(a::AbstractMatrix, b::AbstractMatrix) = kron(a, b)
⊗(a::AbstractVector, b::AbstractVector) = kron(a, b)

function rand_vec(dims...; rng=Random.GLOBAL_RNG)
    return randn(rng, dims...)
end

# Generation of positive-but-not-completely-positive maps. `⊗` and `rand_vec`
# above must be defined before this include so `pncp.jl` can use them.
include("pncp.jl")

# ── Polynomial ↔ Choi matrix conversion ───────────────────────────────────────

function poly2mat(coeffs::AbstractVector, n::Int, m::Int)
    @polyvar X[1:n] Y[1:m]
    xy = X ⊗ Y
    return poly2mat(coeffs ⋅ (xy ⊗ xy), n, m)
end

function poly2mat(poly::AbstractPolynomial, n::Int, m::Int)
    d = n * m
    vars = variables(poly.x)
    X = vars[1:n]
    Y = vars[n+1:end]
    M = zeros(d, d)
    for row in 1:d
        for col in row:d
            i = div(row - 1, m) + 1
            j = div(col - 1, m) + 1
            k = mod(row - 1, m) + 1
            l = mod(col - 1, m) + 1

            mon = X[i] * X[j] * Y[k] * Y[l]
            val = coefficient(poly, mon)

            mult = ((i != j) + 1) * ((k != l) + 1)
            val /= mult

            M[row, col] = val
            M[col, row] = val
        end
    end
    return M
end

# ── Gram-matrix freedom ─────────────────────────────────────────────────────────

"""
    gram_freedom(n, m) -> Vector{Matrix{Float64}}

Basis of the space `L` of symmetric `(n·m)×(n·m)` matrices that vanish on the real
Segre variety: `zᵀ N z ≡ 0` for every product vector `z = x ⊗ y`. The Gram matrix
of a biquadratic form is unique only up to `L`, so `poly2mat`'s output `M₀` and
every representative `Mλ = M₀ + Σ λα·N[α]` describe the *same* polynomial.

The basis encodes the Segre (2×2 minor) relations
`(xᵢyⱼ)(xₖyₗ) − (xᵢyₗ)(xₖyⱼ) ≡ 0`; there are `binomial(n,2)·binomial(m,2)` of them
and `L ≅ ⋀²(ℝⁿ) ⊗ ⋀²(ℝᵐ)`. This is the same construction used inside
`non_sos_form` (the `E(i,j,k,l) - E(i,l,k,j)` span).

Caveat: these vanish on *real* product vectors only. On a complex product vector
`⟨x⊗y| N |x⊗y⟩` is generally nonzero, so adding `L` preserves real block-positivity
but *not* complex block-positivity (the genuine witness condition — see
[`is_block_positive`](@ref)).
"""
function gram_freedom(n::Int, m::Int)
    I_n, I_m = Matrix{Float64}(I, n, n), Matrix{Float64}(I, m, m)
    e(i, j) = kron(I_n[:, i], I_m[:, j])          # basis vector for the monomial xᵢ yⱼ
    basis = Matrix{Float64}[]
    for i in 1:n-1, j in 1:m-1, k in i+1:n, l in j+1:m
        a, b, c, d = e(i, j), e(k, l), e(i, l), e(k, j)
        push!(basis, (a*b' + b*a') - (c*d' + d*c'))
    end
    return basis
end

"""
    is_block_positive(W, n, m; field=:complex, trials=100_000, atol=1e-9, rng) -> Bool

Monte-Carlo check that `W` is block-positive on the `[n, m]` bipartition, i.e.
`⟨x⊗y| W |x⊗y⟩ ≥ -atol` for all product vectors. With `field=:complex` (the genuine
entanglement-witness / positive-map condition) product vectors range over
`ℂⁿ ⊗ ℂᵐ`; with `field=:real` over `ℝⁿ ⊗ ℝᵐ`.

The sampler returns `false` as soon as it finds a violating product vector, so it
can *disprove* block-positivity but never certify it — a certificate needs an
SOS / Positivstellensatz SDP. Note that real block-positivity is invariant under
the [`gram_freedom`](@ref) directions while complex block-positivity is not, so
`field=:real` cannot tell representatives `Mλ` apart.
"""
function is_block_positive(W::AbstractMatrix, n::Int, m::Int;
                           field::Symbol=:complex, trials::Int=100_000,
                           atol::Real=1e-9, rng=Random.GLOBAL_RNG)
    for _ in 1:trials
        z = if field === :complex
            kron(randn(rng, ComplexF64, n), randn(rng, ComplexF64, m))
        else
            kron(randn(rng, n), randn(rng, m))
        end
        real(z' * W * z) < -atol && return false
    end
    return true
end

# ── Map composition ───────────────────────────────────────────────────────────

"""
    ampliation(A, B, n, m; system=2) -> Matrix

Apply the map `Φ_A` (Choi matrix `A`) to one subsystem of `B`, an operator on the
`[n, m] = kron(A_n, B_m)` bipartition (an `(n·m)×(n·m)` matrix).

- `system=2` (default): `(I_n ⊗ Φ_A)(B)`. `Φ_A` acts on the **second**
  (dimension-`m`) subsystem, so `A = J(Φ_A : M_m → M_r)` is `(m·r)×(m·r)` in
  `[m, r]` ordering, `r = size(A,1) ÷ m`, and the result is `(n·r)×(n·r)` in
  `[n, r]` ordering. When `B = J(Φ_B : M_n → M_m)` this is the composition Choi
  matrix `J(Φ_A ∘ Φ_B)` — the form used throughout the PPT² search (square
  self-composition `ampliation(M, M, d, d)`, `n = m = r = d`).

- `system=1`: `(Φ_A ⊗ I_m)(B)`. `Φ_A` acts on the **first** (dimension-`n`)
  subsystem, so `A = J(Φ_A : M_n → M_r)` is `(n·r)×(n·r)`, `r = size(A,1) ÷ n`,
  and the result is `(r·m)×(r·m)` in `[r, m]` ordering. This is the leg dual to a
  block-positive *witness*: for entanglement detection `(Φ_W ⊗ I)(ρ) ⊁ 0` is the
  test consistent with `tr(W·ρ) < 0` (see [`detect_ampliation`](@ref)), whereas
  `system=2` tests the opposite leg and need not agree with the witness.

Rectangular dimensions are supported as long as `Φ_A`'s input dimension matches
the chosen subsystem (`m` for `system=2`, `n` for `system=1`).
"""
function ampliation(A::AbstractMatrix, B::AbstractMatrix, n::Int, m::Int; system::Int=2)
    system == 1 || system == 2 || throw(ArgumentError("system must be 1 or 2, got $system"))
    size(B, 1) == n * m || throw(DimensionMismatch(
        "B must be $(n*m)×$(n*m) for the [$n, $m] bipartition, got $(size(B, 1))×$(size(B, 2))"))

    if system == 2                          # (I_n ⊗ Φ_A)(B): map on subsystem 2 (dim m)
        size(A, 1) % m == 0 || throw(DimensionMismatch(
            "A's dimension $(size(A, 1)) is not a multiple of m=$m (Φ_A's input on subsystem 2)"))
        r = size(A, 1) ÷ m
        # gather B by its second (m) subsystem; natA[(k,l),(a,b)] = A[(a,k),(b,l)]
        Bmat = reshape(PermutedDimsArray(reshape(B, m, n, m, n), (1, 3, 2, 4)), m*m, n*n)
        natA = reshape(PermutedDimsArray(reshape(A, r, m, r, m), (1, 3, 2, 4)), r*r, m*m)
        C = reshape(natA * Bmat, r, r, n, n)            # (k, l, α, β)
        return reshape(PermutedDimsArray(C, (1, 3, 2, 4)), n*r, n*r)
    else                                    # (Φ_A ⊗ I_m)(B): map on subsystem 1 (dim n)
        size(A, 1) % n == 0 || throw(DimensionMismatch(
            "A's dimension $(size(A, 1)) is not a multiple of n=$n (Φ_A's input on subsystem 1)"))
        r = size(A, 1) ÷ n
        # gather B by its first (n) subsystem; natA[(k,l),(a,c)] = A[(a,k),(c,l)]
        Bmat = reshape(PermutedDimsArray(reshape(B, m, n, m, n), (2, 4, 1, 3)), n*n, m*m)
        natA = reshape(PermutedDimsArray(reshape(A, r, n, r, n), (1, 3, 2, 4)), r*r, n*n)
        C = reshape(natA * Bmat, r, r, m, m)            # (k, l, b, d)
        # output in [r, m] ordering: r-leg (k,l) slow, untouched m-leg (b,d) fast
        return reshape(PermutedDimsArray(C, (3, 1, 4, 2)), r*m, r*m)
    end
end

# ── Random states ──────────────────────────────────────────────────────────────

"""
    rand_ppt(n, m; rng, rand_vec) -> Matrix

Random PPT state on the `[n, m] = kron(A_n, B_m)` bipartition (the same ordering
as `rand_sep`, `is_ppt`, and `Ket`). A random PSD matrix whose off-diagonal
m×m blocks (the blocks of the first, dimension-n subsystem) are symmetrised, so
the partial transpose over either subsystem stays PSD; shifted to be PSD if
needed. Pass a custom `rand_vec` to control the entry distribution (e.g. integer
entries).
"""
function rand_ppt(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec, ppt_invariant=false)
    A = rand_vec(n*m, n*m; rng=rng)
    rho = A * A'
    if ppt_invariant
        for i in 1:n, j in i+1:n
            rows = (i - 1) * m + 1:i * m
            cols = (j - 1) * m + 1:j * m
            sym = (rho[rows, cols] + rho[cols, rows]) / 2
            rho[rows, cols] = sym
            rho[cols, rows] = sym
        end
    end
    delta = eigmin(partial_transpose(rho, 2, [n, m]))
    if delta < 0
        return rho - delta * I
    end
    return rho
end

"""Separable state: sum of `n_terms` random product projectors |a⟩⟨a|⊗|b⟩⟨b|."""
function rand_sep(n::Int, m::Int; n_terms::Int=2, rng=Random.GLOBAL_RNG)
    rho = zeros(n * m, n * m)
    for _ in 1:n_terms
        a = randn(rng, n)
        b = randn(rng, m)
        rho += kron(a * a', b * b')
    end
    return rho
end

"""Random PSD matrix of rank `r` (full rank when `r == 0`)."""
function rand_psd(n::Int, m::Int; r::Int=0, rng=Random.GLOBAL_RNG)
    d = n * m
    r = r > 0 ? r : d
    rho = zeros(d, d)
    for _ in 1:r
        psi = randn(rng, d)
        rho += psi * psi'
    end
    return rho
end

# ── Subsystem helpers ──────────────────────────────────────────────────────────

"""Swap operator on C^d ⊗ C^d."""
function swap(d::Int)
    V = zeros(d, d, d, d)
    for i in 1:d, j in 1:d
        V[j, i, i, j] = 1.0
    end
    return reshape(V, d^2, d^2)
end

# Projectors onto the (anti)symmetric subspace of C^d ⊗ C^d. `symmetric_projector`
# is not exported because `Ket` already exports that name; reach it as
# `ppt2.symmetric_projector` when both modules are in scope.
symmetric_projector(d::Int)     = (I(d^2) + swap(d)) / 2
antisymmetric_projector(d::Int) = (I(d^2) - swap(d)) / 2

# ── Smallest-eigenvalue sign checks ───────────────────────────────────────────
#
# Detecting entanglement / non-PPT-ness only needs the SIGN of the smallest
# eigenvalue, not its value. A Cholesky factorisation of `M + tol·I` (LAPACK
# `potrf!`) decides that sign ~6× faster than a full `eigvals` on the 16×16
# matrices of the PPT² search, so it is the primitive behind `is_ppt` and
# `detect_ampliation`. `min_eig` is used only where the value itself is needed.

"""
    has_negative_eig!(M; tol=1e-8) -> Bool

True iff `λ_min(M) < -tol`, via an in-place Cholesky of `M + tol·I` (`potrf!`,
`info ≠ 0` ⟺ not positive-definite). Overwrites `M` (Hermitian/symmetric; only the
upper triangle is read) — for hot loops that rebuild `M` each iteration.
"""
function has_negative_eig!(M::AbstractMatrix; tol::Real=1e-8)
    @inbounds for i in axes(M, 1)
        M[i, i] += tol
    end
    return last(LinearAlgebra.LAPACK.potrf!('U', M)) != 0
end

"""
    has_negative_eig(M; tol=1e-8) -> Bool

`λ_min(M) < -tol`, decided by Cholesky (see [`has_negative_eig!`](@ref)). Works on
a copy, leaving `M` untouched.
"""
has_negative_eig(M::AbstractMatrix; tol::Real=1e-8) = has_negative_eig!(Matrix(M); tol=tol)

"""
    min_eig(M) -> Real

Smallest eigenvalue of the Hermitian/symmetric `M` via the LAPACK selected-range
routine (`eigvals(Hermitian(M), 1:1)`) — only the bottom eigenvalue, no full
spectrum. Use when the value (not just its sign) is needed.
"""
min_eig(M::AbstractMatrix) = eigvals(Hermitian(Matrix(M)), 1:1)[1]

"""True if ρ is PPT: the partial transpose has no eigenvalue below `-tol`. Decided
by the Cholesky sign check [`has_negative_eig`](@ref)."""
function is_ppt(ρ::AbstractMatrix, dA::Int, dB::Int; tol=1e-8)
    PT = partial_transpose(Matrix(ρ), 2, [dA, dB])
    return !has_negative_eig(PT; tol=tol)
end

# ── Entanglement detection ────────────────────────────────────────────────────
#
# Three independent criteria; each returns the raw score, the form/witness
# achieving it, and a `detected` flag. Entangled when trace/ampliation < -tol or
# robustness > tol.

"Linear-witness criterion: min `tr(form·τ)` over `forms`."
# `tr(form·τ)` is real for a real-symmetric witness and Hermitian τ; `real` both
# drops the numerical imaginary part and keeps `findmin` well-defined when τ is
# complex (e.g. the Hermitian PPT states from `min_ppt_witness`).
# `tr(form·τ) = real(dot(form, τ))` for a real-symmetric `form` and Hermitian `τ`
# (tr(AB) = Σ A_ij B_ji, and B_ji = conj(B_ij)), computed WITHOUT materialising the
# 10000 matrix products `forms .* Ref(τ)` did — those allocated ~20 MB per call and,
# in the @threads detection loops, drove constant stop-the-world GC that collapsed
# parallelism. `findmin(f -> …, forms)` scans allocation-free.
detect_trace(τ, forms; tol=1e-8) = let (v, i) = findmin(form -> real(dot(form, τ)), forms)
    (value=v, idx=i, detected=v < -tol)
end

"""
    detect_ampliation(τ, forms, n, m; tol=1e-8, system=1)

Positive-map criterion: scan `forms` for one whose ampliation `(Φ_form ⊗ I)(τ)`
(subsystem `system=1`, the leg dual to the block-positive witness — see
[`ampliation`](@ref)) has an eigenvalue below `-tol`. Detection uses the Cholesky
sign check [`has_negative_eig`](@ref); `value` is the most negative eigenvalue
among the forms that fire (`0.0` if none), computed exactly (via [`min_eig`](@ref))
only for those, and `idx` is the firing form (`0` if none).
"""
function detect_ampliation(τ, forms, n, m; tol=1e-8, system::Int=1)
    best = 0.0
    idx = 0
    for (i, W) in enumerate(forms)
        M = ampliation(W, τ, n, m; system=system)
        has_negative_eig(M; tol=tol) || continue
        v = min_eig(M)
        v < best && (best = v; idx = i)
    end
    return (value=best, idx=idx, detected=idx != 0)
end

"""
    detect_dps(τ, n, m; level=2, tol=1e-8)

Level-`level` DPS robustness from Ket; `witness` is Ket's entanglement witness;
`detected` when robustness `> tol`.

If `τ` is stored complex but is real up to numerical noise — as the witness PPT
states are (`min_ppt_witness` attains its optimum on the real slice, leaving an
imaginary part ~1e-12) — the imaginary part is dropped so `entanglement_robustness`
runs over the real PSD cone instead of the double-dimension complex Hermitian one.
This is the same real-slice trick `test_ppt2.jl` applies on load, here guaranteed at
the solve itself: measured ≈38× faster (2.4 s vs 92 s) with identical robustness.
A genuinely complex `τ` (non-negligible imaginary part) is left untouched.
"""
function detect_dps(τ, n, m; level=2, tol=1e-8)
    M = Matrix(τ)
    if eltype(M) <: Complex && maximum(abs ∘ imag, M; init=0.0) <= 1e-9 * maximum(abs ∘ real, M; init=1.0)
        M = real(M)
    end
    r, w = entanglement_robustness(Hermitian(M), [n, m], level; solver=Mosek.Optimizer)
    return (value=r, witness=w, detected=r > tol)
end

# One named criterion applied to τ.
function _criterion(c::Symbol, τ, forms, n, m, level, tol)
    c === :trace      && return detect_trace(τ, forms; tol=tol)
    c === :ampliation && return detect_ampliation(τ, forms, n, m; tol=tol)
    c === :dps        && return detect_dps(τ, n, m; level=level, tol=tol)
    error("unknown criterion $(c)")
end

"""
    test_ppt2(ρ, σ=ρ; n=4, m=4, compose=true, criteria=(:trace,:ampliation,:dps),
              forms=nothing, level=2, tol=1e-8, mode=:sequential)

Run the detection criteria on a PPT² candidate. With `compose` (the conjecture's
setting) the tested operator is the composite τ = (I⊗Φ_σ)(ρ) =
`ampliation(ρ, σ, n, m)`; `σ` defaults to `ρ` for self-composition. With
`compose=false`, `ρ` itself is tested.

`mode` controls how `criteria` are combined:
- `:sequential` — evaluate in the given order (cheap first), short-circuit on the
  first that fires, and return its evidence `(criterion, value, idx/witness,
  detected)`, or `nothing`. Best for a search loop.
- `:parallel` — evaluate all criteria and return a NamedTuple keyed by criterion
  symbol with each result, plus an overall `detected`. Best for recording every
  score.

`:trace`/`:ampliation` need `forms`.
"""
function test_ppt2(ρ, σ=ρ; n::Int=4, m::Int=4, compose::Bool=true,
                   criteria=(:trace, :ampliation, :dps), forms=nothing,
                   level::Int=2, tol::Float64=1e-8, mode::Symbol=:sequential)
    (:trace in criteria || :ampliation in criteria) && forms === nothing &&
        error("criteria $(criteria) require `forms`")
    τ = compose ? Hermitian(ampliation(ρ, σ, n, m)) : Hermitian(Matrix(ρ))

    if mode === :sequential
        for c in criteria
            d = _criterion(c, τ, forms, n, m, level, tol)
            d.detected && return (criterion=c, d...)
        end
        return nothing
    elseif mode === :parallel
        results = map(c -> _criterion(c, τ, forms, n, m, level, tol), criteria)
        return merge(NamedTuple{criteria}(results),
                     (detected = any(r -> r.detected, results),))
    else
        error("mode must be :sequential or :parallel, got $(mode)")
    end
end

# ── Witness-restricted PPT minimisation ──────────────────────────────────────
#
# The convex dual of `detect_trace`. There the state τ is fixed and we scan a
# finite library of witnesses for the most negative `tr(W·τ)`; here a single
# block-positive witness `W` is fixed and we search the *whole* PPT cone for the
# state it detects most strongly. A negative optimum certifies a PPT entangled
# (bound entangled) state — another way of testing/generating entanglement.

"""
    min_ppt_witness(W, n, m; tol=1e-8, verbose=false, trace_preserving=false)
        -> (value, state, detected)

Minimise `tr(W·ρ)` over operators `ρ` on the `[n, m]` bipartition that are
Hermitian, PSD, and PPT (partial transpose over the second subsystem PSD), via an
SDP. `W` is a fixed `(n·m)×(n·m)` block-positive PnCP witness; `ρ` is the variable.

The scale of `ρ` is fixed one of two ways:

  * `trace_preserving=false` (default): unit trace `tr(ρ) = 1` — `ρ` is a quantum
    *state*, the right object when this SDP is used to *extract a bound entangled
    state* from `W`.
  * `trace_preserving=true`: `tr_2[ρ] = I_n` — `ρ` is the Choi matrix of a
    trace-preserving map (`tr(ρ) = n`). This is the **convex relaxation of the
    see-saw [`min_ppt2_witness`](@ref)**: every composite of two PPT *channels*
    lives in this set, so the optimum here lower-bounds the see-saw optimum on the
    same scale. (Use `n == m` so the Choi matrix matches `W`.)

Every separable state `σ` satisfies `tr(W·σ) ≥ 0`, so a negative optimum
`value < -tol` certifies the minimiser `state` as a PPT *entangled* (bound
entangled) state witnessed by `W`, and sets `detected`. This is dual to
[`detect_trace`](@ref): that fixes the state and ranges the witness over a finite
library; this fixes the witness and ranges the state over the PPT cone.

`state` is returned as a `Hermitian{ComplexF64}`. For a real-symmetric `W` the
optimum is already attained on the real-symmetric slice, so ranging `ρ` over
complex Hermitian operators (as quantum states demand) costs no detection power.
"""
function min_ppt_witness(W::AbstractMatrix, n::Int, m::Int; tol=1e-8, verbose=false,
                         trace_preserving::Bool=false)
    d = n * m
    size(W) == (d, d) || throw(DimensionMismatch(
        "W must be $(d)×$(d) for the [$n, $m] bipartition, got $(size(W))"))

    model = Model(Mosek.Optimizer)
    verbose || set_silent(model)

    @variable(model, ρ[1:d, 1:d] in HermitianPSDCone())                 # ρ ⪰ 0, Hermitian
    if trace_preserving
        @constraint(model, partial_trace(Matrix(ρ), 2, [n, m]) .== Matrix{ComplexF64}(I, n, n))  # trace preserving: tr_2[ρ] = I_n (Choi of a TP map; tr = n)
    else
        @constraint(model, real(tr(ρ)) == 1)                            # unit trace (quantum state)
    end
    @constraint(model, Hermitian(partial_transpose(Matrix(ρ), 2, [n, m])) in HermitianPSDCone())  # PPT
    @objective(model, Min, real(tr(W * ρ)))

    optimize!(model)
    v = objective_value(model)
    return (value=v, state=value.(ρ), detected=v < -tol)
end

# ── Witness-restricted PPT² composition minimisation (see-saw) ────────────────
#
# Like `min_ppt_witness`, but restricts the search to the states that actually
# arise as a *composition* of two PPT maps — the PPT² setting itself — rather than
# the whole PPT cone. The tested operator is `composite = ampliation(ρ1, ρ2)`,
# which is *bilinear* in (ρ1, ρ2): the joint minimisation is a non-convex bilinear
# matrix inequality, not an SDP. We solve it by see-saw (alternating SDPs), the
# same scheme as `min_xy_form`: freeze one factor so the composite is affine in
# the other, minimise that SDP, swap, and repeat from random PPT restarts.

"""
    min_ppt2_witness(W, n, m; restarts=16, max_iter=40, tol=1e-8, verbose=false, rng)
        -> (value, ρ1, ρ2, composite, detected)

See-saw search for a PPT² counterexample detected by a fixed block-positive
witness `W`. Minimise `tr(W · composite)`, where
`composite = ampliation(ρ1, ρ2, n, m)` is the Choi matrix of a composition
`Φ_1 ∘ Φ_2`, over the Choi matrices `ρ1`, `ρ2` of two PPT *channels* on the
`[n, m]` bipartition — Hermitian, PSD, PPT, and trace preserving (`tr_2[ρ] = I_n`).
Trace preservation is closed under composition, so the composite is itself a
PPT channel with `tr(composite) = n`; this fixes the composite's scale and avoids
the trivial degeneracy of an unconstrained-trace composite shrinking to zero.

A composition of PPT maps is itself PPT, so a **negative** optimum
`value < -tol` exhibits a composition `Φ_1 ∘ Φ_2` of PPT maps whose Choi matrix is
PPT *and* entangled (detected by `W`) — i.e. a **counterexample to the PPT²
conjecture**, witnessed by `W`.

`ampliation(ρ1, ρ2, …)` is bilinear in `(ρ1, ρ2)`, so the joint problem is a
non-convex bilinear matrix inequality, not an SDP. It is solved by **see-saw**:
freeze `ρ2` and minimise over `ρ1` (an SDP — the composite is then affine in
`ρ1`), freeze `ρ1` and minimise over `ρ2`, and alternate to convergence from
`restarts` random PPT starts, the same scheme as [`min_xy_form`](@ref). The
returned `value` is therefore only a *local* optimum: a negative one is a genuine
certificate, but a non-negative one does not rule out a counterexample for `W`.
The convex relaxation [`min_ppt_witness`](@ref) — minimising over the *whole* PPT
cone — is a lower bound, so if *it* is `≥ 0` no composition can be negative either.

The composite's own PSD/PPT constraints are dropped because they hold
automatically for a composition of PPT maps; only the two factors are constrained.
"""
function min_ppt2_witness(W::AbstractMatrix, n::Int, m::Int;
                          restarts::Int=16, max_iter::Int=40, tol=1e-8,
                          verbose::Bool=false, rng=Random.GLOBAL_RNG)
    d = n * m
    n == m || throw(ArgumentError(
        "min_ppt2_witness needs n == m (square Choi) so the composite matches W; got n=$n, m=$m"))
    size(W) == (d, d) || throw(DimensionMismatch(
        "W must be $(d)×$(d) for the [$n, $m] bipartition, got $(size(W))"))

    # One persistent PPT-cone model per factor (Φ_1 ⇒ ρ1, Φ_2 ⇒ ρ2). Only the
    # bilinear objective is rebuilt each half-step, with the other factor frozen
    # to its current numeric value, which makes that half an honest SDP.
    function ppt_factor()
        model = Model(Mosek.Optimizer)
        set_silent(model)
        @variable(model, ρ[1:d, 1:d] in HermitianPSDCone())                 # ρ ⪰ 0, Hermitian
        @constraint(model, partial_trace(Matrix(ρ), 2, [n, m]) .== Matrix{ComplexF64}(I, n, n))  # trace preserving: tr_2[ρ] = I_n
        @constraint(model, Hermitian(partial_transpose(Matrix(ρ), 2, [n, m])) in HermitianPSDCone())  # PPT
        return model, ρ
    end
    modelA, ρ1 = ppt_factor()
    modelB, ρ2 = ppt_factor()

    best_val = Inf
    best_ρ1 = Matrix{ComplexF64}(I, d, d) / d
    best_ρ2 = Matrix{ComplexF64}(I, d, d) / d

    for r in 1:restarts
        ρ2cur = rand_ppt(n, m; rng=rng)
        ρ2cur = ρ2cur / tr(ρ2cur)                                           # unit-trace PPT start
        ρ1cur = best_ρ1
        prev = Inf
        for it in 1:max_iter
            # optimise Φ_1 with Φ_2 frozen: composite affine in ρ1
            @objective(modelA, Min, real(tr(W * ampliation(Matrix(ρ1), ρ2cur, n, m))))
            optimize!(modelA)
            ρ1cur = Matrix(value.(ρ1))

            # optimise Φ_2 with Φ_1 frozen
            @objective(modelB, Min, real(tr(W * ampliation(ρ1cur, Matrix(ρ2), n, m))))
            optimize!(modelB)
            ρ2cur = Matrix(value.(ρ2))

            val = objective_value(modelB)
            verbose && println("restart=$r iter=$it value=$val")
            abs(prev - val) < tol && break
            prev = val
        end

        val = real(tr(W * ampliation(ρ1cur, ρ2cur, n, m)))
        if val < best_val
            best_val, best_ρ1, best_ρ2 = val, ρ1cur, ρ2cur
        end
    end

    composite = Hermitian(ampliation(best_ρ1, best_ρ2, n, m))
    return (value=best_val, ρ1=best_ρ1, ρ2=best_ρ2, composite=composite,
            detected=best_val < -tol)
end

# All dataset I/O (readers + the batch-generation engine the scripts drive).
include("io.jl")

end # module ppt2
