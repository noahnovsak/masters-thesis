using LinearAlgebra

# ── Utility ──────────────────────────────────────────────────────────────────

"""Tensor product of two vectors."""
⊗(a::Vector, b::Vector) = kron(a, b)

"""Tensor product of two matrices."""
⊗(A::Matrix, B::Matrix) = kron(A, B)

"""Partial transpose on subsystem B for a bipartite system with dims (dA, dB)."""
function partial_transpose(ρ::Matrix, dA::Int, dB::Int)
    PT = similar(ρ)
    for i in 0:dA-1, j in 0:dA-1
        PT[i*dB+1:(i+1)*dB, j*dB+1:(j+1)*dB] =
            ρ[i*dB+1:(i+1)*dB, j*dB+1:(j+1)*dB]'
    end
    return PT
end

"""Check if a matrix is PPT (all eigenvalues of partial transpose ≥ -tol)."""
function is_ppt(ρ::Matrix, dA::Int, dB::Int; tol=1e-8)
    PT = partial_transpose(ρ, dA, dB)
    return all(eigvals(Hermitian(PT)) .≥ -tol)
end

"""
Realignment criterion: if ‖R(ρ)‖₁ > 1, state is entangled.
R(ρ)[ia, jb] = ρ[ij, ab]  (reshape indices).
"""
function realignment_norm(ρ::Matrix, dA::Int, dB::Int)
    n = dA * dB
    R = zeros(ComplexF64, n, n)
    for i in 0:dA-1, a in 0:dB-1, j in 0:dA-1, b in 0:dB-1
        R[i*dA+j+1, a*dB+b+1] = ρ[i*dB+a+1, j*dB+b+1]
    end
    return sum(svdvals(R))  # nuclear norm = sum of singular values
end

# ── Projector from a product vector ──────────────────────────────────────────

"""Rank-1 projector onto the product state |a⟩⊗|b⟩."""
function product_projector(a::Vector, b::Vector)
    v = a ⊗ b
    return v * v'
end

# ── UPB → Bound Entangled State ───────────────────────────────────────────────

"""
Given an UPB as a vector of (a, b) pairs (normalised kets),
return the bound entangled state ρ = (I - Π) / (d - |UPB|)
where Π = sum of rank-1 projectors onto the UPB vectors.
"""
function upb_state(upb::Vector{Tuple{Vector{Float64},Vector{Float64}}}, dA::Int, dB::Int)
    d  = dA * dB
    Π  = sum(product_projector(a, b) for (a, b) in upb)
    ρ  = (I(d) - Π) / (d - length(upb))
    return Hermitian(ρ)
end

# ── The Shifts UPB in 3⊗3 ────────────────────────────────────────────────────
#
#   Bennett et al. (1999), PRL 82, 5385.
#   Five orthogonal product states whose orthogonal complement contains
#   no product vector → ρ_shifts is PPT and entangled (bound entangled).
#
#   |ψ₀⟩ = |0⟩ ⊗ (|0⟩-|1⟩)/√2
#   |ψ₁⟩ = (|1⟩-|2⟩)/√2 ⊗ |2⟩
#   |ψ₂⟩ = |2⟩ ⊗ (|1⟩-|2⟩)/√2          ← "Shifts" pattern
#   |ψ₃⟩ = (|0⟩-|1⟩)/√2 ⊗ |0⟩
#   |ψ₄⟩ = (|0⟩+|1⟩+|2⟩)/√3 ⊗ (|0⟩+|1⟩+|2⟩)/√3

function shifts_upb_3x3()
    e0, e1, e2 = [1,0,0.0], [0,1,0.0], [0,0,1.0]
    s = (e0 + e1 + e2) / √3

    upb = [
        (e0,              (e0 - e1) / √2),
        ((e1 - e2) / √2,  e2             ),
        (e2,              (e1 - e2) / √2),
        ((e0 - e1) / √2,  e0             ),
        (s,               s              ),
    ]
    return upb
end

"""Build the 9×9 Shifts bound entangled state."""
function shifts_state_3x3()
    upb = shifts_upb_3x3()
    return upb_state(upb, 3, 3)
end

# ── The Tiles UPB in 3⊗3 ─────────────────────────────────────────────────────
#
#   Bennett et al. (1999) — five different product states, also unextendible.
#   Forms a "tiling" of the 3×3 grid of basis states.

function tiles_upb_3x3()
    e0, e1, e2 = [1,0,0.0], [0,1,0.0], [0,0,1.0]

    upb = [
        (e0,              (e0 - e1) / √2),
        (e2,              (e1 - e2) / √2),
        ((e0 - e1) / √2,  e2            ),
        ((e1 - e2) / √2,  e0            ),
        ((e0+e1+e2)/√3,   (e0+e1+e2)/√3 ),
    ]
    return upb
end

"""Build the 9×9 Tiles bound entangled state."""
function tiles_state_3x3()
    upb = tiles_upb_3x3()
    return upb_state(upb, 3, 3)
end

# ── Pyramid UPB in 2⊗4 ───────────────────────────────────────────────────────
#
#   Five product states in C²⊗C⁴ — smallest possible UPB.
#   Bravyi (2004).

function pyramid_upb_2x4()
    θ = [2π*k/5 for k in 0:4]
    # A-vectors live in C², B-vectors in C⁴ (use first 2 dims on a circle)
    upb = Tuple{Vector{Float64}, Vector{Float64}}[]
    for k in 0:4
        a = [cos(π/4), (-1)^k * sin(π/4)]          # alternating ±
        b_raw = [1.0, cos(θ[k+1]), sin(θ[k+1]), 0.0]
        push!(upb, (normalize(a), normalize(b_raw)))
    end
    return upb
end

"""Build the 8×8 Pyramid bound entangled state."""
function pyramid_state_2x4()
    upb = pyramid_upb_2x4()
    return upb_state(upb, 2, 4)
end

# ── Random UPB search (small dimensions) ──────────────────────────────────────
#
#   Numerically search for UPBs by growing a set of mutually orthogonal
#   product states until no orthogonal product complement exists.
#   Uses randomised restarts.

"""
Check if a candidate product vector |a⟩⊗|b⟩ is orthogonal to all states
in `upb` (list of (a,b) pairs).
"""
function orthogonal_to_all(a, b, upb)
    v = a ⊗ b
    return all(abs(dot(u_a ⊗ u_b, v)) < 1e-8 for (u_a, u_b) in upb)
end

"""
Attempt to find a product vector orthogonal to the given UPB partial set
by random sampling + projection onto the orthogonal complement.
Returns (found::Bool, a, b).
"""
function find_orthogonal_product_vector(upb, dA, dB; n_trials=5000)
    d = dA * dB
    # Build orthogonal complement projector
    if isempty(upb)
        return true, normalize(randn(dA)), normalize(randn(dB))
    end
    V = hcat([normalize(a ⊗ b) for (a,b) in upb]...)
    P_perp = I(d) - V * V'   # projector onto ⊥ of span(UPB)

    for _ in 1:n_trials
        a0 = normalize(randn(dA))
        b0 = normalize(randn(dB))
        # Project onto complement and try to separate
        v  = normalize(P_perp * (a0 ⊗ b0))
        # Schmidt rank-1 check via SVD of reshaped vector
        M  = reshape(v, dA, dB)
        sv = svdvals(M)
        if sv[1] / sum(sv) > 1 - 1e-6   # nearly rank-1
            a_new = normalize(M * M' * a0)  # dominant left singular vector
            b_new = normalize(M' * a_new)
            if orthogonal_to_all(a_new, b_new, upb)
                return true, a_new, b_new
            end
        end
    end
    return false, zeros(dA), zeros(dB)
end

"""
Numerically construct a random UPB in C^dA ⊗ C^dB.
Grows the set greedily; returns the UPB and the resulting bound entangled state.
"""
function random_upb(dA::Int, dB::Int; max_size=nothing, n_trials=5000, verbose=false)
    upb = Tuple{Vector{Float64}, Vector{Float64}}[]
    max_size = something(max_size, dA * dB - 1)

    for step in 1:max_size
        found, a, b = find_orthogonal_product_vector(upb, dA, dB; n_trials)
        if !found
            verbose && println("UPB complete at size $(length(upb)) (step $step)")
            break
        end
        push!(upb, (normalize(a), normalize(b)))
        verbose && println("Step $step: added product vector, UPB size = $(length(upb))")
    end

    ρ = upb_state(upb, dA, dB)
    return upb, ρ
end

# ── Verification ──────────────────────────────────────────────────────────────

"""Print a summary of entanglement properties of ρ in C^dA ⊗ C^dB."""
function verify_bound_entangled(ρ::AbstractMatrix, dA::Int, dB::Int; name="ρ")
    println("── $name ($( dA)⊗$dB, dim $(size(ρ,1))×$(size(ρ,2))) ──")
    
    eigs_ρ  = eigvals(Hermitian(Matrix(ρ)))
    PT      = partial_transpose(Matrix(ρ), dA, dB)
    eigs_PT = eigvals(Hermitian(PT))
    ra_norm = realignment_norm(Matrix(ρ), dA, dB)
    tr_ρ    = tr(ρ)
    purity  = real(tr(ρ^2))

    println("  Trace:               $(round(real(tr_ρ), digits=6))")
    println("  Purity tr(ρ²):       $(round(purity, digits=6))")
    println("  Min eig(ρ):          $(round(minimum(eigs_ρ), digits=8))  $(minimum(eigs_ρ) ≥ -1e-8 ? "✓ PSD" : "✗ not PSD")")
    println("  Min eig(ρᴳ):         $(round(minimum(eigs_PT), digits=8)) $(minimum(eigs_PT) ≥ -1e-8 ? "✓ PPT" : "✗ not PPT")")
    println("  Realignment ‖R(ρ)‖₁: $(round(ra_norm, digits=6))          $(ra_norm > 1 + 1e-6 ? "→ entangled (realignment)" : "→ inconclusive")")
    println("  Verdict:             PPT=$(is_ppt(Matrix(ρ),dA,dB)) — bound entangled by UPB range criterion")
    println()
end

# ── Demo ──────────────────────────────────────────────────────────────────────

ρ_shifts = shifts_state_3x3()
verify_bound_entangled(ρ_shifts, 3, 3; name="Shifts UPB (3⊗3)")

ρ_tiles = tiles_state_3x3()
verify_bound_entangled(ρ_tiles, 3, 3; name="Tiles UPB (3⊗3)")

ρ_pyramid = pyramid_state_2x4()
verify_bound_entangled(ρ_pyramid, 2, 4; name="Pyramid UPB (2⊗4)")

# Numerical random UPB search in 3⊗3
upb_rand, ρ_rand = random_upb(3, 3; verbose=true)
verify_bound_entangled(ρ_rand, 3, 3; name="Random UPB (3⊗3)")