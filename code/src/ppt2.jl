module ppt2

using Random
using LinearAlgebra
using JuMP
using MosekTools

using DynamicPolynomials: AbstractPolynomial, @polyvar, variables, coefficient
using SumOfSquares: SOSModel, SOSCone, gram_matrix

import Ket: partial_transpose

export pncp_mat, ampliation, rand_ppt, rand_sep, rand_psd, is_ppt, antisymmetric_projector


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

# ── Map composition ───────────────────────────────────────────────────────────

"""
    ampliation(A, B, n, m) -> Matrix

Choi matrix of a map composition, computed as `(I_n ⊗ Φ_A)(B) = J(Φ_A ∘ Φ_B)`.

`B = J(Φ_B)` is the Choi matrix of `Φ_B : M_n → M_m`, an `(n·m)×(n·m)` matrix in
the `[n, m] = kron(A_n, B_m)` ordering. `Φ_A` is applied to the second
(dimension-`m`) subsystem, so `A = J(Φ_A)` is the Choi matrix of
`Φ_A : M_m → M_r`, an `(m·r)×(m·r)` matrix in `[m, r]` ordering; the output
dimension `r = size(A, 1) ÷ m` is inferred. The result is
`J(Φ_A ∘ Φ_B) : M_n → M_r`, an `(n·r)×(n·r)` matrix in `[n, r]` ordering.

Rectangular dimensions are supported as long as the shared dimension `m`
(`Φ_B`'s output and `Φ_A`'s input) matches — that is what makes the maps
composable. The PPT² search uses the square self-composition
`ampliation(M, M, d, d)` (`n = m = r = d`).
"""
function ampliation(A::AbstractMatrix, B::AbstractMatrix, n::Int, m::Int)
    size(B, 1) == n * m || throw(DimensionMismatch(
        "B must be $(n*m)×$(n*m) for Φ_B: M_$n → M_$m, got $(size(B, 1))×$(size(B, 2))"))
    size(A, 1) % m == 0 || throw(DimensionMismatch(
        "A's dimension $(size(A, 1)) is not a multiple of m=$m (Φ_A's input)"))
    r = size(A, 1) ÷ m

    # gather B into blocks of its second (m) subsystem: Bmat[(a,b), (α,β)]
    Bp = PermutedDimsArray(reshape(B, m, n, m, n), (1, 3, 2, 4))
    Bmat = reshape(Bp, m*m, n*n)

    # natural representation of Φ_A: natA[(k,l), (a,b)] = A[(a,k), (b,l)]
    Ap = PermutedDimsArray(reshape(A, r, m, r, m), (1, 3, 2, 4))
    natA = reshape(Ap, r*r, m*m)

    C = reshape(natA * Bmat, r, r, n, n)            # (k, l, α, β)

    return reshape(PermutedDimsArray(C, (1, 3, 2, 4)), n*r, n*r)
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

"""True if ρ is PPT (partial transpose has no eigenvalue below `-tol`)."""
function is_ppt(ρ::AbstractMatrix, dA::Int, dB::Int; tol=1e-8)
    PT = partial_transpose(Matrix(ρ), 2, [dA, dB])
    return eigmin(Hermitian(PT)) ≥ -tol
end

end # module ppt2
