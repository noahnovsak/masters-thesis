# Generation of positive-but-not-completely-positive (PnCP) maps.
#
# The pipeline produces a quadratic form that is positive on the Segre variety
# (product vectors) but is not a sum of squares — the polynomial counterpart of
# a positive map that is not completely positive. The stages are:
#
#   sample_pncp_form   build one candidate form (f, h)            (unverified)
#   find_pncp_poly     retry until a candidate is positive,       (polynomial)
#                      then rationalize it into a certificate
#   pncp_mat           the same certificate as a Choi matrix
#
# Helpers `⊗` and `rand_vec` are defined in `ppt2.jl` ahead of the include.

"""
    segre_kernel_basis(n::Int, m::Int, Z::Matrix) -> Vector{Matrix}

Compute a basis of the kernels of the generators of the Segre variety, one per
point `Z_i`.

Inputs:
    n, m [Int] - size of the problem.
    Z [Matrix] - random points Z_i on the Segre variety as column vectors.

Outputs:
    W [Vector] - basis of the kernels for each point Z_i.

We loop through the rows and columns. For each 2x2 minor created in this
way, we store the derivative of the minor in a vector to create the
gradient vector of the generators.
"""
function segre_kernel_basis(n::Int, m::Int, Z::Matrix)
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
    non_sos_form(
        n::Int, m::Int, Z::Matrix, W::Vector, h::Matrix;
        rng=Random.GLOBAL_RNG, rand_vec=rand_vec,
    ) -> Vector

Compute a quadratic form that is positive on the Segre variety but is *not* a
sum of squares.

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
function non_sos_form(
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
    solve_sos(n::Int, m::Int, f::Vector, h::Matrix, l=0, fix_gram=false) -> Tuple

Solve the SOS feasibility SDP for the form `δ·f + h⋅h`, maximizing `δ`.

Inputs:
    n, m [Int] - problem dimensions
    f [Vector] - quadratic form
    h [Matrix] - linear forms
    l [Int] - degree of the relaxation (default: 0)
    fix_gram [Bool] - when true, force the smallest Gram eigenvalues to zero and
                      return a 'rationalized' polynomial candidate (default: false)

Output:
    (opt, val) [Tuple] where:
      - opt [Bool] optimization terminates as OPTIMAL and margin exceeds `1e-4`
      - val is either:
          * `value.(poly)` (the optimized polynomial) when `fix_gram == false`
          * `p_hat` (rationalized polynomial) when `fix_gram == true && opt == true`
          * `value.(poly)` and `opt == false` when the rationalized candidate is SOS

Note the return type of `val` depends on `fix_gram`.
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
        return rationalize_certificate(gram_matrix(con), poly, relax, n, m)
    end

    return feasible, value.(poly)
end

"""
    rationalize_certificate(gram, p, r, n, m) -> (is_pncp, poly)

Turn the numerical SOS solution into a rational polynomial certificate.

Zero out the smallest Gram eigenvalues, rebuild the corresponding polynomial,
and double-check that it is still *not* SOS. Returns `(true, p_hat)` when the
rationalized polynomial certifies a positive non-CP map, or `(false, value.(p))`
when it turns out to be SOS.
"""
function rationalize_certificate(gram, p, r, n, m)
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
    sample_pncp_form(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec) -> (f, h)

Sample one candidate for a positive-but-not-completely-positive map.

The result is *unverified*: `f` is built to be non-SOS, but positivity on the
Segre variety still has to be checked (see `find_pncp_poly`).

Inputs:
    n, m [Int] - problem dimensions
    rng [AbstractRNG] - RNG used by `rand_vec`
    rand_vec [Function] - random sampler with signature `rand_vec(dims...; rng=...)`

Outputs:
    (f, h) [Tuple] where:
      - f [Vector] - quadratic form
      - h [Matrix] - linear forms
"""
function sample_pncp_form(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)
    d = n + m - 2
    e = (n - 1) * (m - 1)

    # step 1: random linear form and tensor products
    x = rand_vec(n, e + 1; rng=rng)
    y = rand_vec(m, e + 1; rng=rng)
    z = hcat(kron.(eachcol(x), eachcol(y))...)

    # step 2: construct linear forms h0, ..., hd
    h = nullspace(z') * rand_vec(d, d; rng=rng)
    h0 = nullspace(z[:, 1:end - 1]') * rand_vec(d+1; rng=rng)
    h = hcat(h0, h)

    # step 3: construct (non-SOS) quadratic form f
    w = segre_kernel_basis(n, m, z)
    f = non_sos_form(n, m, z, w, h; rng=rng, rand_vec=rand_vec)

    return f, h
end

"""
    find_pncp_poly(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)

Search for a verified positive-but-not-completely-positive map.

Repeatedly samples candidates with `sample_pncp_form` until one is non-SOS and
positive, then returns the rationalized polynomial certificate.

Inputs:
    n, m [Int] - problem dimensions
    rng [AbstractRNG] - RNG used by `rand_vec`
    rand_vec [Function] - random sampler with signature `rand_vec(dims...; rng=...)`

Outputs:
    poly [Polynomial] - rationalized polynomial certificate when successful
    nothing - if no certificate is found within the retry budget
"""
function find_pncp_poly(n::Int, m::Int; rng=Random.GLOBAL_RNG, rand_vec=rand_vec)
    for attempt in 1:20
        f, h = nothing, nothing
        for construction in 1:20
            f, h = sample_pncp_form(n, m; rng=rng, rand_vec=rand_vec)

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
    poly = find_pncp_poly(n, m; rng=rng, rand_vec=rand_vec)
    if poly === nothing
        return nothing
    end
    return poly2mat(poly, n, m)
end
