module ppt2

using Random
using LinearAlgebra
using JuMP
using MosekTools

using DynamicPolynomials: AbstractPolynomial, @polyvar, variables, coefficient
using SumOfSquares: SOSModel, SOSCone, gram_matrix

export pncp_mat, ampliation, rand_ppt


⊗(a::AbstractMatrix, b::AbstractMatrix) = kron(a, b)
⊗(a::AbstractVector, b::AbstractVector) = kron(a, b)

function rand_vec(dims...; rng=Random.GLOBAL_RNG)
    return randn(rng, dims...)
end

"""
    kernel_basis(n::Int, m::Int, Z::Matrix) -> Vector{Matrix}

Computes a basis of kernels of the generators of the Segre variety.

Inputs:
    n, m [Int] - size of the problem.
    Z [Matrix] - random points Z_i on the Segre variety as column vectors.

Outputs:
    W [Vector] - basis of the kernels for each point Z_i.

We loop through the rows and columns. For each 2x2 minor created in this
way, we store the derivative of the minor in a vector to create the
gradient vector of the generators.
"""
function kernel_basis(n::Int, m::Int, Z::Matrix)
    e = (n - 1) * (m - 1)
    W = Vector{Matrix}(undef, e)

    for col in 1:e
        G = zeros(n*m, n*m*e÷4)
        idx = 1
        for i in 1:n-1, j in 1:m-1
            for k in i+1:n, l in j+1:m
                g = zeros(n * m)
                ij = m * (i - 1) + j
                kl = m * (k - 1) + l
                il = m * (i - 1) + l
                kj = m * (k - 1) + j
                z = Z[:, col]
                g[kl] =  z[ij]
                g[kj] = -z[il]
                g[ij] =  z[kl]
                g[il] = -z[kj]
                G[:, idx] = g
                idx += 1
            end
        end
        W[col] = nullspace(G')
    end

    return W
end

"""
    quadratic_form(
        n::Int, m::Int, Z::Matrix, W::Vector, h::Matrix;
        rng=Random.GLOBAL_RNG, rand_vec=rand_vec,
    ) -> Vector

Computes a quadratic form on the Segre variety, which is not a sum of squares.

Inputs:
    n, m [Int] - problem dimensions
    Z [Matrix] - random points on Segre variety
    W [Vector] - basis of kernels of the generators
    h [Matrix] - computed linear forms h_0, ..., h_d
    rng [AbstractRNG] - RNG used by `rand_vec`
    rand_vec [Function] - random sampler with signature `rand_vec(dims...; rng=...)`

Outputs:
    f [Vector] - computed quadratic form
"""
function quadratic_form(
    n::Int, m::Int, Z::Matrix, W::Vector, h::Matrix;
    rng=Random.GLOBAL_RNG, rand_vec=rand_vec,
)
    e = (n - 1) * (m - 1)
    d = n + m - 1 # d+1
    nm = n * m

    I_n = Matrix(I, n, n)
    I_m = Matrix(I, m, m)
    I_nm = Matrix(I, nm, nm)

    function E(i, j, k, l)
        e_ij = I_n[:, i] ⊗ I_m[:, j]
        e_kl = I_n[:, k] ⊗ I_m[:, l]
        return e_ij ⊗ e_kl + e_kl ⊗ e_ij
    end

    matrices = zeros(e * d + nm * (nm - 1) ÷ 2, nm^2)
    row = 1

    for i in 1:e, j in 1:d
        matrices[row, :] = (Z[:, i] ⊗ W[i][:, j])'
        row += 1
    end

    for i in 1:nm-1, j in i+1:nm
        matrices[row, :] = (I_nm[:, i] ⊗ I_nm[:, j])' - (I_nm[:, j] ⊗ I_nm[:, i])'
        row += 1
    end

    span = zeros(nm^2, nm * e ÷ 4 + d^2)
    col = 1

    for i in 1:n-1, j in 1:m-1
        for k in i+1:n, l in j+1:m
            span[:, col] = E(i, j, k, l) - E(i, l, k, j)
            col += 1
        end
    end

    for i in 1:d, j in 1:d
        span[:, col] = h[:, i] ⊗ h[:, j] + h[:, j] ⊗ h[:, i]
        col += 1
    end

    r = rank(span)
    ker = nullspace(matrices)

    while true
        f = ker * rand_vec(size(ker, 2); rng=rng)
        rank([span f]) > r && return f
    end
end

"""
    solve_sos(n::Int, m::Int, f::Vector, h::Matrix, l=0, zero_g=false) -> Tuple

Checks if the given form is a sum of squares.

Inputs:
    n, m [Int] - problem dimensions
    f [Vector] - quadratic form
    h [Matrix] - linear forms
    l [Int] - degree of the relaxation (default: 0)
    fix_gram [Bool] - Force the smallest Gram eigenvalues to zero and return
                      a 'rationalized' polynomial candidate (default: false)

Output:
    (opt, val) [Tuple] where:
      - opt [Bool] optimization terminates as OPTIMAL and margin exceeds `1e-4`
      - val is either:
          * `del::Float64` when `fix_gram == false`
          * `p_hat` (polynomial) when `fix_gram == true && opt == true`
          * `nothing` when the rationalized candidate is SOS
"""
function solve_sos(n::Int, m::Int, f::Vector, h::Matrix, l=0, fix_gram=false)
    model = SOSModel(Mosek.Optimizer)
    set_silent(model)

    @polyvar X[1:n] Y[1:m]
    @variable(model, δ)

    xy = X ⊗ Y
    f = f ⋅ (xy ⊗ xy)
    h = h'xy

    poly = δ * f + h ⋅ h
    relax = (xy ⋅ xy)^l

    con = @constraint(model, poly * relax in SOSCone())
    @objective(model, Max, δ)
    optimize!(model)

    δ_opt = value(δ)
    feasible = termination_status(model) == OPTIMAL && δ_opt > 1e-4

    if feasible && fix_gram
        return _fix_gram(gram_matrix(con), poly, relax, n, m)
    end

    return feasible, value.(poly)
end

function _fix_gram(gram, p, r, n, m)
    Q = Matrix(gram.Q)
    b = gram.basis.monomials

    vals, vecs = eigen(Q)
    vals[1:(n-1)*(m-1)+1] .= 0
    Q_hat = vecs * Diagonal(vals) * vecs'

    g = b' * Q_hat * b

    A = [coefficient(px * r, gx) for gx in g.x, px in p.x]
    p_hat = A \ g.a ⋅ p.x

    # double check that p_hat is not SOS
    test = SOSModel(Mosek.Optimizer)
    set_silent(test)

    @constraint(test, p_hat in SOSCone())
    @objective(test, Min, 0.0)
    optimize!(test)

    if is_solved_and_feasible(test)
        return false, value.(p)
    end

    return true, p_hat
end

"""
    gen_pncp(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)

Generate a positive but not completely positive map.

Inputs:
    n, m [Int] - problem dimensions
    rng [AbstractRNG] - RNG used by `rand_vec`
    rand_vec [Function] - random sampler with signature `rand_vec(dims...; rng=...)`

Outputs:
    (f, h) [Tuple] where:
      - f [Vector] - quadratic form
      - h [Matrix] - linear forms
"""
function gen_pncp(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)
    d = n + m - 2
    e = (n - 1) * (m - 1)

    # step 1: random linear form and tensor products
    x = rand_vec(n, e + 1; rng=rng)
    y = rand_vec(m, e + 1; rng=rng)
    z = hcat(kron.(eachcol(x), eachcol(y))...)

    # step 2: construct liear forms h0, ..., hd
    h = nullspace(z') * rand_vec(d, d; rng=rng)
    h0 = nullspace(z[:, 1:end - 1]') * rand_vec(d+1; rng=rng)
    h = hcat(h0, h)

    # step 3: construct (non-SOS) quadratic form f
    w = kernel_basis(n, m, z)
    f = quadratic_form(n, m, z, w, h; rng=rng, rand_vec=rand_vec)

    return f, h
end

"""
    pncp_poly(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)

Generate positive maps that are not completely positive.

Inputs:
    n, m [Int] - problem dimensions
    rng [AbstractRNG] - RNG used by `rand_vec`
    rand_vec [Function] - random sampler with signature `rand_vec(dims...; rng=...)`

Outputs:
    poly [Polynomial] - rationalized polynomial certificate when successful
    nothing - if no certificate is found within the retry budget
"""
function pncp_poly(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)
    for attempt in 1:20
        f, h = nothing, nothing
        for construction in 1:20
            f, h = gen_pncp(n, m; rng=rng, rand_vec=rand_vec)

            # verify construction is not SOS
            sos, _ = solve_sos(n, m, f, h)
            if !sos
                break
            end
        end

        # verify construction is positive
        pos, poly = solve_sos(n, m, f, h, 1, true)
        if pos
            return poly
        end
    end

    return nothing
end

function pncp_mat(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)
    poly = pncp_poly(n, m; rng=rng, rand_vec=rand_vec)
    if poly === nothing
        return nothing
    end
    return poly2mat(poly, n, m)
end

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

"""
Computes (I ⊗ A)(B) 
"""
function ampliation(A::AbstractMatrix, B::AbstractMatrix, n::Int, m::Int)
    Ap = PermutedDimsArray(reshape(A, m, n, m, n), (1, 3, 2, 4))
    Bp = PermutedDimsArray(reshape(B, m, n, m, n), (1, 3, 2, 4))

    Ar = reshape(Ap, m*m, n*n)
    Br = reshape(Bp, n*n, n*n)

    Cr = reshape(Ar * Br, m, m, n, n)

    return reshape(PermutedDimsArray(Cr, (1, 3, 2, 4)), n*m, n*m)
end

function rand_ppt(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)
    A = rand_vec(n*m, n*m; rng=rng)
    rho = A * A'
    for i in 1:m, j in i+1:m
        rows = (i - 1) * n + 1:i * n
        cols = (j - 1) * n + 1:j * n
        sym = (rho[rows, cols] + rho[cols, rows]) / 2
        rho[rows, cols] = sym
        rho[cols, rows] = sym
    end
    delta = eigmin(rho)
    if delta < 0
        return rho - delta * I
    end
    return rho
end

end # module ppt2
