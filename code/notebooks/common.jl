# Notebook-only helpers: form quality-control and the known bound-entangled
# example. `include("common.jl")` once at the top of a notebook.
#

using Random
using LinearAlgebra
using JuMP
using MosekTools
using Ket
using JLD2
using ProgressMeter

using ppt2

# For integer-entry PPT states, pass a custom sampler to rand_ppt, e.g.
#     rand_ppt(n, m; rand_vec = (d...; rng) -> float(rand(rng, -1:1, d...)))

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
