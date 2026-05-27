# Shared helpers for the experiment notebooks.
#
# Include once at the top of a notebook with `include("common.jl")`.
# Method-specific generators (UPB search, antisymmetric-subspace SDP) live in
# their own notebooks; this file holds the pieces that were being copy-pasted
# between several of them: the PPT2 detection driver, form quality-control,
# the random-state zoo, and a few small projectors.

using Random
using LinearAlgebra
using JuMP
using MosekTools
using Ket
using JLD2
using ProgressMeter

using ppt2   # gen_pncp, pncp_mat, solve_sos, ampliation, rand_ppt, ...

# ── Loading precomputed PNCP forms ────────────────────────────────────────────

"""
    load_forms(path) -> Vector{Matrix}

Concatenate every batch stored in a `pncp_forms_*.jld2` file into one vector.
"""
function load_forms(path::AbstractString)
    jldopen(path, "r") do file
        vcat([file[k] for k in keys(file) if !startswith(k, "meta/")]...)
    end
end

# ── Random-state zoo ──────────────────────────────────────────────────────────
#
# `rand_ppt` itself lives in the `ppt2` library and is the canonical generator.
# To sample PPT states with integer entries (as in the old `gen_ppt` notebook),
# pass a custom sampler, e.g.
#
#     rand_ppt(n, m; rand_vec = (d...; rng) -> float(rand(rng, -1:1, d...)))

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

"""Haar-random unitary of size `n` via QR with phase correction."""
function random_unitary(n::Int; rng=Random.GLOBAL_RNG)
    X = (randn(rng, n, n) .+ im * randn(rng, n, n)) ./ sqrt(2)
    F = qr(X)
    Q = Matrix(F.Q)
    phases = [iszero(F.R[i, i]) ? 1.0 + 0im : F.R[i, i] / abs(F.R[i, i]) for i in 1:n]
    return Q * Diagonal(phases)
end

# ── A known 3⊗3 bound-entangled example ───────────────────────────────────────

"""
    example_bound_entangled(a) -> 9×9 Matrix

The one-parameter Horodecki 3⊗3 bound entangled state (`0 ≤ a ≤ 1`). PPT for all
`a`, entangled for `0 < a < 1`. Useful as a known-positive sanity check.
"""
function example_bound_entangled(a::Float64)
    @assert 0.0 <= a <= 1.0
    aii = (1 + a) / 2
    aij = sqrt(1 - a^2) / 2
    return 1 / (8 * a + 1) * [
          a  0.0  0.0  0.0    a  0.0  0.0  0.0    a
        0.0    a  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0    a  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0    a  0.0  0.0  0.0  0.0  0.0
          a  0.0  0.0  0.0    a  0.0  0.0  0.0    a
        0.0  0.0  0.0  0.0  0.0    a  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  aii  0.0  aij
        0.0  0.0  0.0  0.0  0.0  0.0  0.0    a  0.0
          a  0.0  0.0  0.0    a  0.0  aij  0.0  aii
    ]
end

"""Anti-diagonal ±1 witness operator on C^d (detects the example above)."""
function B(d::Int)
    b = zeros(d, d)
    for k in 1:d
        b[k, end - k + 1] = (-1)^k
    end
    return b
end

# ── Symmetric / antisymmetric subspace projectors ─────────────────────────────

"""Swap operator on C^d ⊗ C^d."""
function swap(d::Int)
    V = zeros(d, d, d, d)
    for i in 1:d, j in 1:d
        V[j, i, i, j] = 1.0
    end
    return reshape(V, d^2, d^2)
end

symmetric_projector(d::Int)     = (I(d^2) + swap(d)) / 2
antisymmetric_projector(d::Int) = (I(d^2) - swap(d)) / 2

"""True if ρ is PPT (partial transpose has no eigenvalue below `-tol`)."""
function is_ppt(ρ::AbstractMatrix, dA::Int, dB::Int; tol=1e-8)
    PT = partial_transpose(Matrix(ρ), 2, [dA, dB])
    return eigmin(Hermitian(PT)) ≥ -tol
end

# ── PNCP form quality control ─────────────────────────────────────────────────
#
# A generated form `f` should define a *positive* map: ⟨x⊗y| f |x⊗y⟩ ≥ 0 for all
# product vectors. These two routines look for a product vector that violates it.

"""
    is_pos(form, n, m; attempts, tol) -> (min_value, x, y)

Random product-vector search for a negative value of ⟨x⊗y| form |x⊗y⟩.
Returns the smallest value seen and the witnessing (x, y) once it drops below
`-tol` (else `(min_value, nothing, nothing)`).
"""
function is_pos(form, n::Int, m::Int; attempts=100_000, tol=1e-6, rng=Random.GLOBAL_RNG)
    mi = Inf
    for _ in 1:attempts
        x = randn(rng, n)
        y = randn(rng, m)
        xy = kron(x, y)
        mi = min(mi, xy' * form * xy)
        mi < -tol && return mi, x, y
    end
    return mi, nothing, nothing
end

"""
    min_xy_form(form, n, m; restarts, max_iter, tol) -> (min_value, x, y)

Minimize ⟨x⊗y| form |x⊗y⟩ by an alternating ("see-saw") SDP: optimize the X
factor with Y fixed, then Y with X fixed, repeating to convergence from several
random restarts. A negative minimum certifies the form is *not* positive.
"""
function min_xy_form(
    form, n::Int, m::Int; restarts::Int=16, max_iter::Int=40,
    tol::Float64=1e-8, verbose::Bool=false, rng=Random.GLOBAL_RNG,
)
    Q = Matrix(Hermitian(form))

    best_val = Inf
    best_x = zeros(n)
    best_y = zeros(m)

    modelX = Model(Mosek.Optimizer); set_silent(modelX)
    modelY = Model(Mosek.Optimizer); set_silent(modelY)
    @variable(modelX, X[1:n, 1:n], PSD)
    @constraint(modelX, sum(X[i, i] for i in 1:n) == 1.0)
    @variable(modelY, Y[1:m, 1:m], PSD)
    @constraint(modelY, sum(Y[j, j] for j in 1:m) == 1.0)

    for r in 1:restarts
        y0 = normalize(randn(rng, m))
        Ycur = y0 * y0'
        Xcur = nothing

        prev = Inf
        for it in 1:max_iter
            A = zeros(n, n)
            for i in 1:n, k in 1:n
                A[i, k] = sum(Q[(i - 1) * m + j, (k - 1) * m + l] * Ycur[j, l]
                              for j in 1:m, l in 1:m)
            end
            @objective(modelX, Min, sum(A[i, k] * X[i, k] for i in 1:n, k in 1:n))
            optimize!(modelX)
            Xcur = Symmetric(value.(X))

            B = zeros(m, m)
            for j in 1:m, l in 1:m
                B[j, l] = sum(Q[(i - 1) * m + j, (k - 1) * m + l] * Xcur[i, k]
                              for i in 1:n, k in 1:n)
            end
            @objective(modelY, Min, sum(B[j, l] * Y[j, l] for j in 1:m, l in 1:m))
            optimize!(modelY)
            Ycur = Symmetric(value.(Y))

            val = sum(Q[(i - 1) * m + j, (k - 1) * m + l] * Xcur[i, k] * Ycur[j, l]
                      for i in 1:n, j in 1:m, k in 1:n, l in 1:m)
            verbose && println("restart=$r iter=$it value=$val")
            abs(prev - val) < tol && break
            prev = val
        end

        ex = eigen(Hermitian(Xcur))
        ey = eigen(Hermitian(Ycur))
        x = ex.vectors[:, argmax(ex.values)]
        y = ey.vectors[:, argmax(ey.values)]
        val_pure = real(kron(x, y)' * Q * kron(x, y))
        if val_pure < best_val
            best_val, best_x, best_y = val_pure, x, y
        end
    end

    return best_val, best_x, best_y
end

# ── PPT2 detection ────────────────────────────────────────────────────────────
#
# Three independent ways to certify that a state τ is entangled:
#   :trace      — minimal ⟨form, τ⟩ over the PNCP forms (a linear witness)
#   :ampliation — minimal eigenvalue of (I⊗form)(τ) over the forms
#   :dps        — DPS robustness from Ket (level-2 hierarchy)
# A negative trace / ampliation value, or a positive robustness, means entangled.

detect_trace(τ, forms; tol=1e-8) = let (v, i) = findmin(tr.(forms .* Ref(τ)))
    (value=v, idx=i, detected=v < -tol)
end

detect_ampliation(τ, forms, n, m; tol=1e-8) =
    let (v, i) = findmin(minimum.(real.(eigvals.(ampliation.(forms, Ref(τ), n, m)))))
        (value=v, idx=i, detected=v < -tol)
    end

detect_dps(τ, n, m; tol=1e-8) = let (r, w) = entanglement_robustness(Hermitian(Matrix(τ)), [n, m], 2; solver=Mosek.Optimizer)
    (value=r, witness=w, detected=r > tol)
end

"""
    test_ppt2(generator; n=4, m=4, n_trials=1000, compose=true,
              criteria=(:trace, :ampliation, :dps), forms=nothing, tol=1e-8)

Search for a counterexample to the PPT² conjecture.

`generator(; rng)` must return a PPT state ρ ∈ C^n ⊗ C^m. When `compose` is true
(the conjecture's setting) the candidate tested is the composite
τ = (I⊗Φ_ρ)(ρ) built by `ampliation(ρ, ρ, n, m)`; otherwise ρ itself is tested.
Each requested criterion is checked; `:trace`/`:ampliation` need `forms`.

Returns `(state, evidence)` on the first detection, else `(nothing, nothing)`.
"""
function test_ppt2(
    generator; n::Int=4, m::Int=4, n_trials::Int=1000, compose::Bool=true,
    criteria=(:trace, :ampliation, :dps), forms=nothing, tol::Float64=1e-8,
    rng=Random.GLOBAL_RNG,
)
    if (:trace in criteria || :ampliation in criteria) && forms === nothing
        error("criteria $(criteria) require `forms` (load with load_forms)")
    end

    @showprogress for _ in 1:n_trials
        ρ = generator(; rng=rng)
        τ = compose ? Hermitian(ampliation(ρ, ρ, n, m)) : Hermitian(Matrix(ρ))

        if :trace in criteria
            d = detect_trace(τ, forms; tol=tol)
            d.detected && return ρ, (criterion=:trace, d...)
        end
        if :ampliation in criteria
            d = detect_ampliation(τ, forms, n, m; tol=tol)
            d.detected && return ρ, (criterion=:ampliation, d...)
        end
        if :dps in criteria
            d = detect_dps(τ, n, m; tol=tol)
            d.detected && return ρ, (criterion=:dps, d...)
        end
    end
    return nothing, nothing
end
