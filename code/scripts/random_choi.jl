using LinearAlgebra, Random

function random_unitary(n::Int)
    X = (randn(n, n) .+ im * randn(n, n)) ./ sqrt(2)
    F = qr(X)
    Q = Matrix(F.Q)
    R = F.R
    phases = similar(diag(R))
    for i in eachindex(phases)
        r = R[i, i]
        phases[i] = iszero(r) ? 1.0 : r / abs(r)
    end
    return Q * Diagonal(phases)
end

function random_entangled_choi(dim_in::Int=2, dim_out::Int=2)
    d_total = dim_in * dim_out

    # random unitary on the total space
    U = random_unitary(d_total)

    # construct maximally entangled vector (match Python's row-major flattening)
    phi_plus = zeros(ComplexF64, d_total)
    for i in 1:min(dim_in, dim_out)
        idx = (i - 1) * dim_out + i
        phi_plus[idx] = 1
    end
    phi_plus ./= sqrt(dim_in)

    # apply unitary
    psi = U * phi_plus

    # pure state projector
    rho = psi * psi'

    # interpret as Choi and normalize trace to dim_in
    choi = rho .* (dim_in / tr(rho))

    return choi
end

function partial_transpose(rho::AbstractMatrix{T}, dim_in::Int, dim_out::Int) where T
    d = dim_in * dim_out
    rho_pt = similar(rho)
    for i in 1:dim_in, a in 1:dim_out, k in 1:dim_in, b in 1:dim_out
        orig_row = (i - 1) * dim_out + a
        orig_col = (k - 1) * dim_out + b
        new_row = (i - 1) * dim_in + k
        new_col = (a - 1) * dim_out + b
        rho_pt[new_row, new_col] = rho[orig_row, orig_col]
    end
    return rho_pt
end

# Example usage when run as a script
if abspath(PROGRAM_FILE) == @__FILE__
    choi_matrix = random_entangled_choi(2, 2)
    println(choi_matrix)

    eigvals = eigen(Hermitian(choi_matrix)).values
    println("Eigenvalues: ", eigvals)
    println("Is PSD? ", all(eigvals .>= -1e-12))

    rho_pt = partial_transpose(choi_matrix, 2, 2)
    eigvals_pt = eigen(Hermitian(rho_pt)).values
    println("Partial transpose eigenvalues: ", eigvals_pt)
    println("Is entangled? ", any(eigvals_pt .< -1e-12))
end
