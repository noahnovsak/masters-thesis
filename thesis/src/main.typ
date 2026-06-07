#import "conf.typ": *

#show: style-algorithm
#show: checklist
#show: thm-rules

#show: conf.with(
  title_en: "A Software Approach to the\nPPT2 Conjecture",
  title_sl: "Programski pristop k domnevi PPT2",
  author: "Noah Novšak",
  mentor: "doc. dr. Aljaž Zalar",
  cosupervisor: "dr. Igor Klep",
  keywords_en: "PPT2 conjecture, quantum entanglement, positive maps, semidefinite programming, bound entanglement",
  keywords_sl: "domneva PPT2, kvantna prepletenost, pozitivne preslikave, semidefinitno programiranje, vezana prepletenost",
  code_url: "https://github.com/noahnovsak/masters-thesis",
  abstract_en: [This thesis explores...],
  abstract_sl: [V tem delu raziskujemo...],
  extended_abstract_sl: [Daljši slovenski povzetek vsebine...],
)

// fixes "layout did not converge" warning
// #let margin-note = margin-note.with(dy: 0pt)

// TEMPORARY: red placeholder for a result number still to be filled in.
// Grep for `nbr(` to find all pending numbers; remove this helper once none remain.
#let nbr(x) = text(fill: red)[‹#x›]

= TODO <todo>

*Writing*
- [ ] Write the abstracts: English, Slovene, and the extended Slovene summary.
- [ ] Polish the Introduction: spell out the implications of the conjecture holding versus failing, and sharpen the statement of current knowledge ($n <= 3$ proven, $n >= 4$ open).
- [ ] Expand the Conclusion into a genuine discussion: the Gram-representation freedom, the inertness of the symmetric (partial transpose invariant) family, the limitations, and the outlook.

*Results* (pending)
- [ ] Finalize the PPT2 search results: final candidates and witnesses, hardware, wall-clock, and outcome.
- [ ] Fill in the performance numbers in the Performance section: Julia generation timings and the MATLAB-prototype comparison.
- [ ] Remove the temporary `nbr()` red-highlight helper once all numbers are final.

*References*
- [ ] Complete and verify all citations: resolve the remaining reference placeholders and margin notes.


= Introduction <intro>

Quantum entanglement is a central nonclassical resource in quantum information theory. It supports communication and cryptographic tasks and provides a structural lens for understanding quantum channels @Horodecki_2009.

This thesis studies the PPT2 conjecture: whether the composition of two PPT maps is always entanglement breaking @Christandl_2019. The conjecture has a clear operational interpretation, since entanglement-breaking maps destroy all bipartite entanglement with any reference system.

The stakes are structural. Were the conjecture to hold, composing two PPT maps would always destroy entanglement: the PPT cone, though strictly larger than the set of entanglement-breaking maps, would collapse into it after a single self-composition -- a strong rigidity statement about how far a PPT map can sit from being entanglement breaking. Were it to fail, the counterexample would be a pair of PPT maps whose composition is itself PPT yet entangled, a bound entangled channel assembled from ordinary PPT ingredients, bearing on the structure of quantum channels and on the distillability questions that turn on the PPT property @Horodecki_2009. Either way the answer sharpens the picture of the boundary between positive and completely positive maps.

Known progress motivates a computation-first strategy. The conjecture is proven for $n <= 3$ and for Choi-type maps in all dimensions @Chen_2019 @Singh_2022, but it remains open for $n >= 4$, where the geometry of positive maps grows markedly more intricate and direct analytic classification has so far resisted. The smallest open case, $4 times 4$, is the natural target for a computational attack, and the one we pursue.

The core challenge is twofold. First, separability testing is hard in the regimes where the conjecture remains open. Second, numerical certificates near feasibility boundaries are fragile. Therefore, this thesis focuses on methodology: how to construct, compose, test, and validate candidate objects in a reproducible way.

The main objectives are:
1. Design a Julia workflow for candidate generation, map composition, and entanglement checks.
2. Integrate DPS-based semidefinite relaxations as a baseline witness route.
3. Integrate an SOS-based PNCP witness-construction route with post-solver validation.

The rest of the thesis is organized as follows. Chapter 2 develops the mathematical background: notation, matrix spaces, linear maps, entanglement, the criteria used to detect it, and the semidefinite-programming machinery, including the polynomial route to witnesses. Chapter 3 presents the computational methods and the implementation. Chapter 4 reports the experimental results. Chapter 5 discusses their significance, the limitations, and directions for future work.

= Theoretical Background

This chapter establishes the mathematical framework underpinning the PPT2 conjecture and its computational study. We fix notation and define the key objects: linear maps, their positivity properties, quantum states, and entanglement criteria, then develop the semidefinite programming tools used in the implementation.

== Notation <notation>

We write $CC^n$ and $RR^n$ for the complex and real coordinate spaces. $M_n (FF)$ is the algebra of $n times n$ matrices over the field $FF in {RR, CC}$, abbreviated $M_n$ when the field is clear from context. A linear map between matrix algebras is written $Phi: M_n -> M_m$. We write $E_(i j) in M_n$ for the matrix unit with a $1$ in position $(i,j)$ and zeros elsewhere.

In some parts of this work we also use _Dirac notation_, as it is the standard convention in quantum mechanics. With that in mind, let $|psi chevron.r$ represent the column vector $psi in CC^n$, and $chevron.l psi|$ its adjoint row vector $psi^*$. Now $chevron.l psi | phi chevron.r = psi^* phi in CC$ is the inner product and $|psi chevron.r chevron.l psi| = psi psi^* in M_n$ is the rank-one projector onto $psi$; in particular $E_(i j) = |i chevron.r chevron.l j|$. A bipartite _product vector_ is a simple tensor $|psi chevron.r times.o |phi chevron.r$, also written $|psi times.o phi chevron.r$ or $|psi phi chevron.r$.

== Matrix spaces

We equip $M_n$ with involution $*$, that is, conjugate transposition ($A^* = overline(A^T)$) over $CC$ and transposition ($A^* = A^T$) over $RR$. We denote the subspace of Hermitian matrices ($A = A^*$) by $H_n subset.eq M_n (CC)$, and the subspaces of symmetric ($A = A^T$) and skew-symmetric ($A = -A^T$) matrices by $S_n, K_n subset.eq M_n (RR)$ respectively, so that $M_n (RR) = S_n plus.o K_n$. Finally, $M_n^+ subset.eq H_n$ is the cone of positive semidefinite matrices.

#definition(name: "Positive semidefinite and positive definite")[
  A matrix $A in H_n$ is _positive semidefinite_ (PSD), written $A succ.eq 0$, if $x^* A x >= 0$ for all $x in CC^n$, or equivalently if all eigenvalues of $A$ are non-negative. It is _positive definite_ (PD), written $A succ 0$, if $x^* A x > 0$ for all $x in CC^n without {0}$, equivalently if all eigenvalues are strictly positive.
]

The _Hilbert-Schmidt inner product_ on $M_n$ is defined as $chevron.l A, B chevron.r := tr(B^* A),$ making $M_n$ an inner product space. In the case of $H_n$ or $S_n$, it simplifies to $tr(B A)$.

We work with bipartite systems on the tensor product vector space $CC^m times.o CC^n$. We equip the operator space $M_m times.o M_n tilde.eq M_(m n)$ with multiplication, i.e. $(A times.o B)(C times.o D) = A C times.o B D$, for $A, C in M_m$ and $B, D in M_n$. We use the _partial trace_ $tr_k$, a generalization of the trace operation, to _trace out_ the $k$-th factor, for example $tr_1 (A times.o B) = tr(A) B$.

#definition(name: "Partial transpose")[
  For a product matrix $A times.o B in M_m times.o M_n$, the _partial transpose_ with respect to subsystem $B$ is defined as
  $ (A times.o B)^(Gamma_B) := A times.o B^T, $
  and extended by linearity to a general $rho in M_m times.o M_n$.
]

For Hermitian $rho$, the partial transpose with respect to either subsystem yields the same result: $rho^(Gamma_A) = rho^(Gamma_B)$, so we abbreviate to $rho^Gamma$ without loss of generality.

== Linear maps

Given linear maps $Phi_1: M_n -> M_m$ and $Phi_2: M_p -> M_q$, their _tensor product_ $Phi_1 times.o Phi_2: M_n times.o M_p -> M_m times.o M_q$ is the linear map fixed on simple tensors by
$ (Phi_1 times.o Phi_2)(A times.o B) = Phi_1 (A) times.o Phi_2 (B) $
and extended by linearity. We write $I_k$ for the identity map on $M_k$, and we say $I times.o Phi$ is the _ampliation_ of $Phi$.

#definition(name: "Positive, k-positive, and completely positive maps")[
  Let $Phi: M_n -> M_m$ be a linear map.
  - $Phi$ is _positive_ (P) if $A succ.eq 0 => Phi(A) succ.eq 0$.
  - $Phi$ is _$k$-positive_ if $I_k times.o Phi: M_k times.o M_n -> M_k times.o M_m$ is positive.
  - $Phi$ is _completely positive_ (CP) if it is $k$-positive for every $k in NN$ @Chen_2019.
]

Clearly every completely positive map is also positive. However, there are many maps that are positive but not completely positive (PNCP). They will be playing a central role in this work.

#definition(name: "Block positivity")[
  An operator $W in H_(m n)$ is _block positive_ if $ chevron.l x times.o y|W|x times.o y chevron.r >= 0 $ for all product vectors $x in CC^m$ and $y in CC^n$.
]

#definition(name: "Choi-Jamiolkowski isomorphism")[
  Let $Phi: M_n -> M_m$ be a linear map, then the ampliation of $Phi$ on the maximally entangled state $|Omega chevron.r = sum_(i=1)^n |i chevron.r times.o |i chevron.r$ is its _Choi matrix_:
  $ C_Phi = (I times.o Phi)(|Omega chevron.r chevron.l Omega|) = sum_(i,j=1)^n |i chevron.r chevron.l j| times.o Phi(|i chevron.r chevron.l j|) in M_n times.o M_m. $
  The assignment $Phi arrow.bar C_Phi$ is a linear isomorphism between maps $Phi: M_n -> M_m$ and matrices in $M_n times.o M_m$. The map can be recovered from its Choi matrix by $Phi(rho) = tr_1 [(rho^T times.o I_m) C_Phi]$.
]

#theorem[
  Under the Choi-Jamiołkowski isomorphism, $Phi$ is completely positive if and only if $C_Phi succ.eq 0$ @Choi_1975, and $Phi$ is trace-preserving if and only if $tr_2 [C_Phi] = I_n$.
]

Block positivity is the operator counterpart of map positivity: a map $Phi$ is positive if and only if its Choi matrix $C_Phi$ is block positive, exactly as it is completely positive if and only if $C_Phi succ.eq 0$.

The partial transpose is itself the ampliation of the transposition map $T$: $rho^Gamma = (I_m times.o T)(rho)$. Writing this out for the $2 times 2$ case we get its Choi matrix:
$ C_T = sum_(i,j=1)^2 E_(i j) times.o E_(i j)^T = mat(1, 0, 0, 0; 0, 0, 1, 0; 0, 1, 0, 0; 0, 0, 0, 1), $
with eigenvalues $plus.minus 1$. Since $C_T succ.eq.not 0$, the transpose is _not_ completely positive even though it is positive; it is the prototypical PNCP map.

A word to the wise: the Choi matrix of the composition $Phi compose Psi$ is not the product $C_Phi C_Psi$. Instead, it is obtained as the ampliation $I_n times.o Phi$ of $C_Psi$:
$ C_(Phi compose Psi) = (I_n times.o Phi)(C_Psi). $ <map-comp>
Expanding in the standard basis yields the index formula
$ (C_(Phi compose Psi))_(i p, j q) = sum_(k,l) (C_Psi)_(i k, j l) (C_Phi)_(k p, l q). $
This formula is used directly in the implementation to compose two candidate PPT maps from their Choi matrices.

== Entanglement <entanglement>

A _quantum state_ on $CC^n$ is a density matrix $rho in M_n^+$ with $tr(rho) = 1$. A _quantum channel_ is a completely positive trace-preserving (CPTP) map $Phi: M_n -> M_m$.

#definition(name: "Separability and entanglement")[
  A bipartite state $rho in M_m times.o M_n$ is _separable_ if
  $ rho = sum_i p_i rho_i^A times.o rho_i^B, quad p_i >= 0, quad sum_i p_i = 1, quad rho_i^A in M_m^+, quad rho_i^B in M_n^+. $
  Otherwise $rho$ is _entangled_.
]

#definition(name: "Entanglement-breaking map")[
  A map $Phi: M_n -> M_m$ is _entanglement breaking_ (EB) if $(I_k times.o Phi)(rho)$ is separable for every $k in NN$ and every state $rho in M_k times.o M_n$ @Horodecki_2009.
] <eb-def>

Testing $k = 1$ alone would only check that $Phi$ preserves separability of bipartite states; the ampliation over all $k$ ensures $Phi$ destroys entanglement with any external reference system. Equivalently, $Phi$ is EB if and only if its Choi matrix $C_Phi$ is a separable state @Horodecki_2009.

#definition(name: "PPT state and PPT map")[
  A state $rho in M_m times.o M_n$ is _PPT_ if $rho succ.eq 0$ and $rho^Gamma succ.eq 0$. A map $Phi: M_n -> M_m$ is _PPT_ if its Choi matrix is PPT:
  $ C_Phi succ.eq 0 quad "and" quad C_Phi^Gamma succ.eq 0. $
]<ppt-def>

Since $C_Phi succ.eq 0$ characterizes CP maps, a PPT map is a CP map whose Choi matrix has non-negative partial transpose. The composition of two PPT maps is trivially PPT.

=== PPT2 Conjecture

With PPT maps (@ppt-def) and entanglement-breaking maps (@eb-def) in hand, we can now state the central object of study.

#theorem(name: "PPT2 Conjecture")[
  If $Phi_1$ and $Phi_2$ are PPT maps, then $Phi_1 compose Phi_2$ is entanglement breaking @Christandl_2019.
]<ppt2>

The conjecture is proven for $n = 2$ (since all PPT states in $M_2 times.o M_n$ are separable @Horodecki_1996), for $n = 3$ @Chen_2019, and for Choi-type maps in all dimensions @Singh_2022. It remains open for $n >= 4$ in the general case. The rest of this chapter develops the formalism needed to study the conjecture computationally.

== Testing entanglement <qse-section>

Deciding separability is NP-hard in general @Gharibian_2009. Since an exact test is computationally intractable, one relies on necessary conditions: a state that fails such a condition must be entangled. However, satisfying all known criteria does not guarantee separability. For a comprehensive survey see @Guhne_2009.

=== The PPT test

If a state $rho$ is separable, i.e. $rho = rho_A times.o rho_B$, then it is also PPT, since $rho^Gamma = rho_A times.o rho_B^T succ.eq 0$. Furthermore, the transposition map $T$ is positive but not completely positive, applying $I times.o T$ to an entangled state can produce a non-PSD result. This means that if $rho^Gamma succ.eq.not 0$, then it must be entangled. The condition is necessary but not sufficient: bound entangled (PPT entangled) states exist for systems with $m n > 6$ @Horodecki_2009.

=== Positive maps and entanglement witnesses

The set of separable states $"SEP"$ is convex and closed. By the Hahn-Banach separation theorem; any point outside a closed convex set is separated from it by a hyperplane, for every entangled $rho in.not "SEP"$ there exists a Hermitian operator $W$ (a hyperplane in $H_(m n)$) with $tr[W rho] < 0$ and $tr[W sigma] >= 0$ for all $sigma in "SEP"$ @Horodecki_2009.

#definition(name: "Entanglement witness")[
  A Hermitian operator $W in H_(m n)$ is an _entanglement witness_ if $tr[W sigma] >= 0$ for every separable state $sigma in M_m times.o M_n$ and $tr[W rho] < 0$ for some entangled state $rho in M_m times.o M_n$.
]

Geometrically, each witness $W$ cuts off a half-space that contains $rho$ but not $"SEP"$; the intersection of all such half-spaces recovers $"SEP"$ exactly.

Witnesses are one-sided: $tr[W rho] < 0$ certifies entanglement, but $tr[W rho] >= 0$ does not certify separability. An operator is a (valid) entanglement witness precisely when it is block positive but not PSD.

#definition(name: "Decomposable witness")[
  A witness $W$ is _decomposable_ if $W = P + Q^Gamma$ for some $P, Q succ.eq 0$; otherwise it is _non-decomposable_ @Lewenstein_2000.
]<decomposable-def>

A known result is that decomposable witnesses cannot detect PPT entangled states. For any PPT $sigma$, we have $ tr[(P + Q^Gamma) sigma] = tr[P sigma] + tr[Q sigma^Gamma] >= 0. $ This means only non-decomposable witnesses are useful for our purposes.

#figure(
  image("figures/entanglement_witness_hyperplanes.svg", width: 80%),
  caption: [Geometric picture (schematic) of entanglement witnesses as separating hyperplanes over the nested convex sets of separable, PPT and entangled states. Inclusions are strict, with bound-entangled states populating the outside of the PPT region (green strip). The entanglement witness $W_1$ (solid) is _finer_ than $W_2$ (dashed), it detects every state $W_2$ does, in addition to the states $W_2$ misses (shaded band).]
)

Under the Choi-Jamiolkowski isomorphism, every entanglement witness $W = C_Phi$ corresponds to a PNCP map $Phi$, and vice versa; block positivity of $W$ is exactly positivity of $Phi$. This gives two ways to test a state. The scalar test asks whether $tr[C_Phi rho] < 0$, i.e. whether $rho$ lies on the negative side of the single separating hyperplane defined by $C_Phi$. The map test asks the more demanding question whether the ampliation $(I_k times.o Phi)(rho)$ has _any_ negative eigenvalue, $(I_k times.o Phi)(rho) succ.eq.not 0$. The scalar test inspects only one expectation value and so can miss a negative eigenvalue in a direction other than the one $C_Phi$ singles out; the map test inspects the whole spectrum. Concretely, $(I_k times.o Phi)(rho) succ.eq.not 0$ holds exactly when some vector $psi$ satisfies
$ chevron.l psi|(I times.o Phi)(rho)|psi chevron.r = tr[(I times.o Phi^*)(|psi chevron.r chevron.l psi|) rho] < 0. $
The map therefore encodes an entire family of scalar witnesses $W_psi = (I times.o Phi^*)(|psi chevron.r chevron.l psi|)$, one per output vector $|psi chevron.r$; the Choi-matrix test $tr[C_Phi rho]$ is the single member obtained from the maximally entangled $|psi chevron.r$, whereas the eigenvalue check optimizes over all $|psi chevron.r$ at once. Hence
$ lambda_min ((I times.o Phi^*)(rho)) <= tr[C_Phi rho], $
so the map detects every state the scalar test does, and in general strictly more @Horodecki_2009.

The PPT test is the special case of this positive-map test in which $Phi = T$ is the transposition map: $rho^Gamma succ.eq.not 0$ is exactly $(I times.o T)(rho) succ.eq.not 0$. The transpose is, however, a _decomposable_ witness, so it sees only NPT entanglement; detecting PPT entanglement requires stronger maps.

=== Further separability criteria

Several other necessary conditions for separability are known; each can detect entanglement the PPT test misses, but none is sufficient.

*Reduction criterion.* If $rho in M_m times.o M_n$ is separable, then
$ rho_A times.o I_n - rho succ.eq 0 quad "and" quad I_m times.o rho_B - rho succ.eq 0, $
where $rho_A = tr_B [rho]$ and $rho_B = tr_A [rho]$ are the _reduced states_ (obtained by tracing out one subsystem) @Horodecki_1999.

*Range criterion.* If $rho in M_m times.o M_n$ is separable, there exist product vectors ${|psi_i times.o phi_i chevron.r}$ spanning the range of $rho$ such that ${|psi_i times.o overline(phi_i) chevron.r}$ spans the range of $rho^Gamma$ @Horodecki_1997. The range criterion can detect certain PPT entangled states, but fails when $rho$ is full rank (e.g., under noise), since then any set of vectors spans its range trivially.

*Majorization criterion.* If $rho$ is separable, then $lambda(rho) prec.eq lambda(rho_A)$ and $lambda(rho) prec.eq lambda(rho_B)$, where $lambda(dot)$ denotes the non-increasingly ordered eigenvalue vector @Nielsen_2001. Note, here $a prec.eq b$ has nothing to do with positive definiteness; it denotes _majorization_: $sum_(i=1)^k a_i <= sum_(i=1)^k b_i$ for all $k$, with equality at $k = n$. This criterion follows from the reduction criterion @Nielsen_2001 and therefore shares its limitations.

*CCNR / realignment criterion.* Define the _realignment_ of $rho$ as the matrix $R(rho)$ with entries $(R(rho))_(i (m+k), j (n+l)) = rho_(i (n+j), k (n+l))$. If $rho$ is separable, then $||R(rho)||_1 <= 1$, where $||dot||_1$ is the trace norm @Chen_2003. The CCNR criterion is independent of the PPT criterion and can detect some PPT entangled states that the other criteria miss.

Beyond criterion-based tests, there are _algorithmic_ approaches that reformulate separability as a convex optimization problem, i.e. the DPS hierarchy (section @dps-section). There are also criteria based on covariance matrices and Bell inequalities that are less directly applicable in our setting, for a full survey see @Guhne_2009.

== Semidefinite Programming

A _semidefinite program_ (SDP) is a convex optimization problem in which a linear objective is minimized subject to a linear matrix inequality. The standard primal form is

$ "minimize" &quad c^T bold(x) \
  "subject to" &quad F(bold(x)) := F_0 + sum_(i=1)^n x_i F_i succ.eq 0, $ <sdp-primal>

where $bold(x) in RR^n$ is the optimization variable, $c in RR^n$ is the cost vector, and $F_0, F_1, ..., F_n in H_d$ are fixed Hermitian matrices. The associated _dual SDP_ is

$ "maximize" &quad -tr[F_0 Z] \
  "subject to" &quad Z succ.eq 0, \
               &quad tr[F_i Z] = c_i, quad i = 1, ..., n. $ <sdp-dual>

The _Slater condition_ for @sdp-primal requires the existence of a strictly feasible point: some $bold(x)$ with $F(bold(x)) succ 0$. When it holds for both primal and dual, _strong duality_ holds: primal and dual optima coincide. When $c = 0$, @sdp-primal is a _feasibility problem_. If infeasible, a dual feasible $Z succ.eq 0$ with $tr[F_i Z] = 0$ and $tr[F_0 Z] > 0$ certifies infeasibility.

=== Interior point methods

SDPs are solved in practice by _interior point methods_ (IPMs), which follow a smooth trajectory through the strict interior of the feasible region. A standard IPM applied to a $d times d$ SDP with $n$ variables reaches $epsilon$-accuracy in $O(sqrt(d) log(1 slash epsilon))$ iterations, each requiring $O(n^2 d^2 + n d^3)$ operations @Vandenberghe_1996. This thesis uses MOSEK @MOSEK, a state-of-the-art IPM solver for semidefinite programs.

A property of IPMs important for our use: solutions lie in the _strict interior_ of the feasible region $F(bold(x)) succ 0$. Since all solutions are floating-point approximations, this provides the numerical slack around the boundary that enables rationalization as a post solver method of certifying solutions.

=== Sum-of-squares relaxation

Semidefinite programming also underlies a general relaxation of polynomial non-negativity. Let $RR[bold(x), bold(y)]$ be the ring of real polynomials in $bold(x) in RR^n$, $bold(y) in RR^m$. A polynomial $p$ is _non-negative_ if $p(bold(x), bold(y)) >= 0$ for all real inputs; it is a _sum of squares_ (SOS) if $p = sum_i q_i^2$ for polynomials $q_i$. Every SOS polynomial is non-negative, but not every non-negative polynomial is SOS.

#definition(name: "Gram matrix representation")[
  Let $bold(z) = (z_1, ..., z_N) in RR^N$ and let $bold(v)(bold(z))$ represent the $binom(N + d - 1, d)$ monomials of degree $d$ in $bold(z)$. A homogeneous polynomial $p$ of degree $2d$ is SOS if and only if there exists a symmetric PSD matrix $G succ.eq 0$ (a _Gram matrix_ of $p$), such that
  $ p(bold(z)) = bold(v)(bold(z))^T G bold(v)(bold(z)). $ <gram-rep>
]

Testing whether $p$ is SOS therefore reduces to a semidefinite feasibility problem: find $G succ.eq 0$ satisfying the linear constraints that equate coefficients of @gram-rep with those of $p$. If no such $G$ exists, $p$ is not SOS.

Testing the nonnegativity of a polynomial $p$ is an NP-hard problem @Luo_1998. The standard way around this is to use SOS relaxations which are computationally tractable. The idea is to multiply $p$ by a fixed SOS multiplier and test the product for SOS. #margin-note[reference, explanation of Positivstellensatz, formula for the multiplier/retrieving polynomial from the multiplier]

=== Entanglement testing as an SDP

We can rephrase entanglement testing as maximizing $tr[M rho]$, where $rho in "SEP"$. This is a convex optimization problem, but the set $"SEP"$ is not semidefinite representable. Therefore, we find tractable outer approximations of $"SEP"$, to get (a sequence of) necessary conditions for separability. Each approximation is defined by semidefinite constraints, so membership can be tested by an SDP feasibility problem. If $rho$ fails the test at some level, we can conclude it must be entangled.

=== The DPS hierarchy <dps-section>

The Doherty-Parrilo-Spedalieri (DPS) hierarchy @Doherty_2004 provides such a sequence of SDP relaxations of separability, with the added benefit of being _complete_: for every entangled state there exists a finite level $k$ at which the test detects it. In principle, running the hierarchy to convergence solves separability exactly; in practice, only the first few levels are computationally feasible.

*Symmetric extensions.* For a separable state $rho = sum_i p_i rho_i^A times.o rho_i^B$, define the $k$-fold _symmetric extension_
$ rho_k = sum_i p_i rho_i^A times.o (rho_i^B)^(times.o k) in H(CC^m times.o (CC^n)^(times.o k)). $
This extension satisfies three properties, each with a physical interpretation:
$ (I_m times.o Pi_k) rho_k (I_m times.o Pi_k) = rho_k, $ <ext-sym>
$ rho_k^(Gamma_s) succ.eq 0 quad forall s = 1,...,k, $ <ext-ppt>
$ tr_(B_2 ... B_k)[rho_k] = rho. $ <ext-marg>
Here $Pi_k$ projects onto the bosonic (fully symmetric) subspace of $(CC^n)^(times.o k)$. Condition @ext-sym states that $rho_k$ is invariant under permutations of the $k$ copies of $B$. Condition @ext-ppt states that all partial transposes of $rho_k$ are PSD. Condition @ext-marg states that tracing out the extra $k-1$ copies of $B$ recovers $rho$.

#definition(name: "DPS set at level k")[
  $"DPS"_n^k$ is the set of states $rho in H(CC^m times.o CC^n)$ for which there exists a symmetric extension $rho_k in H(CC^m times.o (CC^n)^(times.o k))$ satisfying @ext-sym, @ext-ppt, and @ext-marg.
]

Each $"DPS"_n^k$ is defined by semidefinite constraints, so membership is testable in polynomial time for fixed $k$. The hierarchy satisfies @Doherty_2004:
1. $"SEP"_n subset.eq "DPS"_n^k$ and $"DPS"_n^(k+1) subset.eq "DPS"_n^k$ for all $k >= 1$.
2. $"DPS"_n^1$ is equivalent to the PPT criterion.
3. Asymptotic completeness: $inter.big_(k >= 1) "DPS"_n^k = "SEP"_n$.

*Feasibility and witness extraction.* Testing $rho in "DPS"_n^k$ amounts to searching for $rho_k$ satisfying @ext-sym, @ext-ppt, and @ext-marg: this is an SDP feasibility problem. When it is infeasible (no valid extension exists), the dual variable $Z succ.eq 0$ yields an entanglement witness. Concretely, the SDP dual to the level-$k$ test has a feasible $Z$ whenever $rho$ is entangled at that level. The operator
$ W = tr_(B_2 ... B_k) [Z] $
is an entanglement witness for $rho$ @Doherty_2004, i.e. $tr[W rho] < 0$ verifies entanglement directly.

*Improvements and limitations.* Several enhancements to the basic hierarchy are known. Harrow, Natarajan, and Wu @Harrow_2017 add first-order optimality (KKT) conditions to the DPS SDP, achieving _finite convergence_: infeasibility is certified at a finite level rather than only asymptotically. The KKT conditions are linear constraints on the Lagrange multipliers of the original SDP; they increase the variable count substantially and make clean witness extraction from the dual more involved. Specialized hierarchies have been developed for states with symmetry. For diagonal unitary invariant states, Britz and Laurent @Britz_2025 give a drastically smaller SDP at each level. For Werner states and isotropic states, explicit separability conditions are known. However, none of these apply to our search: the PPT2 conjecture is already proven for the relevant symmetric families @Singh_2022, so our search must use the general DPS hierarchy.

=== Maps as polynomials <sec:maps-as-polynomials>

An alternative approach to detecting entanglement exploits the polynomial representation of linear maps. Each linear map $Phi: M_n -> M_m$ corresponds to a biquadratic polynomial:
$ p_Phi (bold(x), bold(y)) := bold(y)^T Phi(bold(x) bold(x)^T) bold(y). $
The fundamental correspondence @Klep_2017 between the map and polynomial representation is:
- $Phi$ is _positive_ if and only if $p_Phi$ is non-negative on $RR^n times RR^m$.
- $Phi$ is _completely positive_ if and only if $p_Phi$ is SOS.
A PNCP map corresponds exactly to a non-negative non-SOS polynomial. Each such a polynomial therefore gives rise to an entanglement witness, and the SOS test of the previous section becomes a test of complete positivity. This is the route we use to manufacture entanglement witnesses directly.

==== KMSZ construction for PNCP maps <kmsz>

Now that we have established the value of non-negative non-SOS polynomials for our purposes, we present a construction @Klep_2017 that produces such polynomials, and therefore PNCP maps, as follows:
1. Sample random points $x^((1)), ..., x^((t)) in RR^n$ and $y^((1)), ..., y^((t)) in RR^m$.
2. Form bilinear $h_j (bold(x), bold(y)) = chevron.l x^((j)), bold(x) chevron.r dot chevron.l y^((j)), bold(y) chevron.r$, each a product of two linear forms, so $sum_j h_j^2$ is SOS.
3. Find $f in.not "span"{h_1, ..., h_t}$ such that $f$ is not SOS.
4. Find $delta > 0$ small enough that $F_delta := delta f + sum_j h_j^2 >= 0$.

Clearly $F_delta$ is by construction _not_ SOS, so all that is left is to verify the non-negativity. Steps 1--3 use only linear algebra, and step 4 involves solving an SDP. However, we cannot directly represent $F_delta$ with semidefinite constraints. Instead, we relax the condition and search for a Gram matrix $G succ.eq 0$ for $F_delta dot S$, where $S = sum_(i,j) x_i^2 y_j^2$ is a fixed SOS multiplier. This bihomogeneous multiplier preserves the separate degrees in $bold(x)$ and $bold(y)$, matching the biquadratic structure of $F_delta$. If feasible, then $F_delta dot S$ is SOS, which certifies that $F_delta$ is non-negative @Klep_2017.

==== Real maps and their complexification <complexification>

The correspondence above is stated for real symmetric inputs: $p_Phi (bold(x), bold(y)) = bold(y)^T Phi(bold(x) bold(x)^T) bold(y)$ probes $Phi$ only on rank-one _real_ symmetric matrices $bold(x) bold(x)^T$, and "non-negative" means $p_Phi >= 0$ on $RR^n times RR^m$. Accordingly, the KMSZ construction produces a _real_ linear map
$ Phi: S_n (RR) -> S_m (RR). $
Quantum states, however, are complex Hermitian operators, so to use $Phi$ as an entanglement witness it must be extended to a complex map. Following @Klep_2017, we take the _complexified trivial extension_
$ Gamma_CC := (Phi plus.o 0)_CC : M_n (CC) -> M_m (CC), $ <complexified-trivial-extension>
defined in two steps. First, using $M_n (RR) = S_n (RR) plus.o K_n (RR)$, extend $Phi$ to all real matrices by acting as zero on the skew-symmetric part:
$ (Phi plus.o 0)(S + K) = Phi(S), quad S in S_n (RR), space K in K_n (RR). $
Then extend $dagger$-linearly to complex matrices,
$ Gamma_CC (A + i B) = (Phi plus.o 0)(A) + i (Phi plus.o 0)(B), quad A, B in M_n (RR), $
the unique extension satisfying $Gamma_CC (X^*) = Gamma_CC (X)^*$. Klep et al. prove that $Gamma_CC$ is positive (respectively PNCP) whenever $Phi$ is @Klep_2017.

Two observations make this extension harmless to compute with. First, for a Hermitian input $rho = A + i B$ (with $A in S_n (RR)$ symmetric and $B in K_n (RR)$ skew-symmetric, since $rho = rho^*$) only the symmetric part survives, $Gamma_CC (rho) = Phi(A)$. Second, and equivalently, the Choi matrix of $Gamma_CC$ is exactly the real symmetric matrix $C_Phi$ returned by the construction, now read as a complex Hermitian operator; no entries change. A by-product of annihilating the skew part is that $C_Phi$ is _partial-transpose invariant_, $C_Phi^Gamma = C_Phi$, a property we return to in section @witness-findings.

==== Indecomposability of the generated maps <indecomposability>

For complex maps the appropriate polynomial object is _bi-Hermitian_: with $z in CC^n$, $w in CC^m$ and their conjugates, the polynomial
$ p_(Gamma_CC) (z, w) = chevron.l w|Gamma_CC (|z chevron.r chevron.l z|)|w chevron.r $
is real-valued, and it is non-negative on all product vectors if and only if $C_(Gamma_CC)$ is block-positive, i.e. an entanglement witness. For such polynomials there are _two_ inequivalent notions of sum-of-squares @Fang_2020:
- _complex SOS_ (CSOS): $p = sum_i |q_i (z, w)|^2$ with each $q_i$ holomorphic (a polynomial in $z, w$, not $overline(z), overline(w)$);
- _real SOS_ (RSOS): $p = sum_i g_i^2$ with each $g_i$ a Hermitian (real-valued) polynomial in $z, overline(z), w, overline(w)$.
CSOS is a strrictly stronger condition than RSOS. Under the map correspondence @Fang_2020, a CSOS certificate corresponds exactly to a completely positive map and an RSOS certificate to a _decomposable_ witness. Failure of any RSOS decomposition is therefore the polynomial signature of indecomposability, and by @decomposable-def such a witness can detect PPT entangled states.

The two representations are linked by a structural identity. Writing $z = a + i b$ and $w = c + i d$ with $a, b in RR^n$ and $c, d in RR^m$, a direct expansion gives
$ p_(Gamma_CC) (a + i b, c + i d) = p_Phi (a, c) + p_Phi (b, c) + p_Phi (a, d) + p_Phi (b, d). $ <complexification-identity>
The complexified polynomial is thus completely controlled by the values of the real form $p_Phi$ on real inputs, which can be used to prove _indecomposability_ @Masse_2026.

#theorem(name: "Indecomposability")[
  If $p_Phi$ is not SOS, then $p_(Gamma_CC)$ is not RSOS, and hence the witness $C_(Gamma_CC)$ is non-decomposable.
]<indecomp-thm>

Since $p_Phi$ is bihomogenous it vanishes whenever either argument is zero, so restricting to the real slice, i.e. setting $b, d = 0$ we find @complexification-identity collapses to $p_(Gamma_CC) (a, c) = p_Phi (a, c)$ and $p_(Gamma_CC) in "RSOS"$ would imply $p_Phi in "SOS".$

Because the KMSZ form $p_Phi$ is non-SOS by construction, every generated witness is indecomposable, and therefore capable of detecting some PPT entangled state. This settles, at least in principle, our concerns about the efficacy of using such maps for our purposes.

Knowing the generated witnesses _can_ detect PPT entanglement is not our only concern; we also need to know evaluate their detection power. This is determined by where they lie in the space of entanglement witnesses.

We say a witness $W$ is _tight_ @Chruscinski_2014 if there exists some product vector $x times.o y$ for which $chevron.l x times.o y|W|x times.o y chevron.r = 0$, geometrically this means the supporting hyperplane genuinely touches $"SEP"$. Further, $W$ is _optimal_ @Lewenstein_2000 if there is no $P succ.eq 0$ such that $W - P$ remains block-positive; i.e. _tilting_ the hyperplane toward the $"SEP"$.
By construction, the generated witnesses are tight: the algorithm prescribes zeros at the sampled Segre points, so the hyperplane defined by $W$ touches $"SEP"$ at those points. Optimality however, is a stronger condition, and as of yet undetermined for these witnesses @Masse_2026. A practical consequence is that the witnesses lie close to the decomposable cone, and while they are robust in detecting NPT entanglement, they detect PPT entanglement only weakly.

==== Non-uniqueness of the matrix representation <gram-freedom>

The polynomial $p_Phi$ does not determine the matrix $C_Phi$ uniquely. Writing $p_Phi$ as a quadratic form in the product monomials $bold(z) = bold(x) times.o bold(y)$,
$ p_Phi (bold(x), bold(y)) = bold(z)^T M bold(z), $
we see any symmetric $M$ that reproduces the coefficients of $p_Phi$ is an admissible representation. Two such matrices $M$, $M'$ give the same polynomial precisely when $bold(z)^T (M - M') bold(z) = 0$. The matrices with this property form a linear space $L$, spanned by the $2 times 2$ minor (_Segre_) relations
$ (x_i y_k)(x_j y_l) - (x_i y_l)(x_j y_k) = 0, quad i < j, space k < l. $
There is one independent relation for each choice of two rows $i < j$ and two columns $k < l$, so $dim L = binom(n, 2) binom(m, 2)$. The admissible representations therefore form an affine space $M_0 + L$.

The important distinction is that the relations defining $L$ vanish on _real_ product vectors, not on complex ones, so members of $M_0 + L$ are not necessarily valid witnesses or even block-positive over $CC$. The canonical choice is the symmetric, partial-transpose-invariant representative @complexified-trivial-extension, where zeroing the skew-symmetric part guarantees positivity. Whether the skew-symmetric directions can be used while preserving positivity over $CC$ remains open.

==== Rationalization <rationalization>

It is possible to sample over $QQ$ to produce a rational polynomial, but the final step in our construction relies on solving an SDP, which inherently works with floating-point arithmetic. For an exact certificate it is therefore more prudent to rationalize only the final result.

Let $G$ be a numerical solution with $mu = min("eig"(G_0)) > 0$ and residual $epsilon = max_i |chevron.l A_i, G_0 chevron.r - b_i|$, then for $mu > epsilon$, a rational feasible $hat(G)$ can be obtained @Peyrl_2008,@Cafuta_2015 in two steps.
1. Compute a rational approximation $hat(G)$ of the gram matrix $G$ satisfying $||hat(G) - G_0||^2 + epsilon^2 < mu^2$.
2. Project $hat(G)$ back to the affine subspace given by the constraints.

With a rational gram matrix $hat(G)$ in hand, the coefficients of $hat(F_delta)$ can be isolated using exact computations.

= Methods <methods>

This chapter defines the computational workflow used in the current implementation. The goal is reproducibility with mathematical traceability. Each computational step corresponds to a defined object or operation from the preceding chapter.

A note on normalization is in order. A quantum state is conventionally a PSD operator with unit trace, and a quantum channel is trace preserving. The tests in this thesis are, however, positivity and sign based. Detecting entanglement amounts to checking the sign of a partial transpose, a witness expectation, or an SOS margin, none of which depend on the magnitude of the trace. We therefore work throughout with unnormalized operators, and normalize only when necessary. This avoids redundant rescaling in the inner loops without affecting any conclusion.

Two design constraints shape the workflow:
1. _exact_ separability testing is intractable in the dimensions where the conjecture is open (@dps-section), partial testing is very much doable, but every complete test is a one-sided relaxation;
2. SDP solutions are floating-point and fragile near feasibility boundaries, so positive detections are re-validated exactly (@rationalization).

== Generating a witness library <witness-library>

Stage 1 builds the witness library. We run the KMSZ construction to produce $N_W$ PNCP maps, and rationalize each certificate so that every stored witness is exact, guarding against the false positives that floating-point SOS tests invite near the boundary.

The rationalization, developed in section @rationalization, proceeds as follows:
1. Extract the Gram matrix $G$ from the SDP constraints, $ (delta f + h^2)(sum_(i,j) x_i^2 y_j^2)^l = (bold(x)bold(y))^T G bold(x) bold(y) $.
2. Eigendecompose $G$ exploiting the known $e$-dimensional nullspace @Klep_2017; set the first $e = (n-1)(m-1)$ eigenvalues to zero, and rationalize the remaining coefficients to obtain $hat(G)$.
3. Reformulating the constraints in terms of $hat(G)$, isolate the coefficients of the polynomial $hat(p) = delta f + h^2$ from the relaxation terms $(sum_(i,j) x_i^2 y_j^2)^l$.
4. Additionally we solve a final feasibility SDP to double check we have not made $hat(p)$ into a SOS. 

Upon completion we are confident in the validity of the certificate, so we can choose to represent it in floating point or as a rational number, as long as we are sure to use high enough precision. Typically, solvers work with a tolerance of $approx 10^(-8)$.

== Generating candidates

We generate PPT states directly as the Choi matrices of PPT maps. A short sampling routine produces a random PPT state; an entanglement filter then keeps only the candidates that can plausibly yield a counterexample. A second, witness-guided generator targets bound entangled states directly.

=== Random PPT states

#algorithm-figure(
  "Generating a random PPT state",
  {
    Line[Sample a real matrix $R in M_(n m)$ with i.i.d. standard-normal entries and form the PSD matrix $rho = R R^T$.]
    Line[Compute the smallest eigenvalue of the partial transpose, $lambda = lambda_min (rho^Gamma)$. If $lambda < 0$, set $rho <- rho - lambda I$.]
  }
) <ppt-gen>

The construction is correct without any further work. By definition $rho = R R^T succ.eq 0$, and when $lambda < 0$ the shift adds the non-negative multiple $|lambda| I$ to $rho$, so positivity is preserved. The same shift acts on the partial transpose as $rho^Gamma - lambda I$, because $I^Gamma = I$, which raises every eigenvalue of $rho^Gamma$ by $|lambda|$ and hence makes $rho^Gamma succ.eq 0$. The resulting $rho$ is therefore PPT. Sampling elemnt-wise from the normal distribution produces a representative collection of states @Zyczkowski_2011.

The off-diagonal blocks may optionally be symmetrized before the shift, for each block pair $(i, j)$ with $i < j$, both $rho_(i j)$ and $rho_(j i)$ are replaced by their average $(rho_(i j) + rho_(j i)) slash 2$, which forces $rho = rho^Gamma$, invariant under partial transpose (IPT). This symmetrization is not essential. We only check if it makes the computation any simpler or faster, since a single positivity shift then certifies both $rho succ.eq 0$ and $rho^Gamma succ.eq 0$ at once. The Choi matrix of a map is only defined up to the basis convention, so restricting to such IPT representatives is a permissible convenience rather than a requirement. The motivation is that our precomputed witnesses are themselves generated in this IPT shape, so matching the candidates to the same form may make them easier to detect; we therefore run the search in two versions, with and without symmetry.

A random PPT state may be separable or bound entangled, and only the entangled ones are useful here. If either composed map is entanglement breaking, equivalently, has a separable Choi matrix, then the composition is automatically entanglement breaking and cannot violate the conjecture. We therefore discard separable candidates and keep only states for which we can certify entanglement, assembling a pool of genuine bound entangled states whose compositions are worth testing.

=== Bound entangled states from a witness <gen-witness-ppt>

The random sampler has no control over entanglement, so many of its draws are discarded. A second generator instead targets entanglement directly, using the witness library to manufacture states it is guaranteed to detect. Fix a witness $W$ and minimize its expectation over the whole PPT cone,
$ "minimize" quad & tr[W rho] \ "subject to" quad & rho succ.eq 0, \ &rho^Gamma succ.eq 0, \ & tr[rho] = 1. $ <min-ppt-witness>
This is a single SDP. Since every separable state gives $tr[W rho] >= 0$, a negative optimum exhibits a PPT _entangled_ state that $W$ detects, one bound entangled candidate per witness that admits one. The optimum also measures the witness's detection strength, and is a ready-made hard instance for the following composition tests.

=== Alternative bound-entangled constructions

The first method, i.e. random sampling is fast, but it gives no control over entanglement. So we must filter states out, which necessitates solving an SDP. Not a great option for scalable searches. The second method is more targeted, but still suffers from the same problem, though the SDP is much smaller in this case. There are several constructions in the literature that produce bound entangled states by design, however we opted not to use them in the search, for various reasons noted below. They are however, the natural starting points for further investigation.

- *Unextendible product bases (UPB)* are a set of mutually orthogonal product vectors spanning a proper subspace whose orthogonal complement contains no product vector. The normalized projector onto that complement is PPT and entangled @Bennett_1999. UPBs give explicit, low-rank bound entangled states, but the construction is the opposite of generic: it relies on specific, hand-picked bases that cannot be sampled at random, and in some dimensions may not exist at all, so it cannot drive a randomized search.
- *Antisymmetric-subspace states* are a simple class of PPT entangled states from the projectors onto the symmetric and antisymmetric subspaces of two identical systems, generalizing the Werner states @Sindici_2018. The construction is explicit, but certifying entanglement of the resulting states still reduces to an SDP, so it does not scale better than our pipeline.
- *Symmetric random induced states* are a more recently studied class of bound entangled states @Louvet_2025. Sampling a random pure state from the symmetric subspace $cal(H)_S^(N + N_A + 1)$ and tracing out $N_A$ subsystems, leaving a mixed state on $cal(H)_S^(N + 1)$ produces, with high probability, a bound entangled state. The catch is dimensionality: for $N = 4$ bound entanglement is most likely at $N_A = 12$, which demands an enormous ancillary space and correspondingly heavy computation. Even then the probability of entanglement stays below $0.5$, comparable to the hit rate of our own construction. It is nonetheless a promising state-of-the-art source of candidates and a useful point of comparison.

== Testing candidates

Given two PPT candidates $Phi_1$, $Phi_2$ from the pool, we form the Choi matrix of their composition using the ampliation operation @map-comp. Then we test the composite for entanglement: if any composition is ever found to be entangled, the PPT2 conjecture is violated. Because composition is not commutative, the search ranges over _every ordered pair_ of pool states, self-pairs included.

=== Composition search with the witness library and DPS

Each composite is screened by three complementary detection tests:
- *Scalar witness test.* For each witness $W$ in the precomputed library, $tr[W C] < 0$ certifies that $C$ is entangled.
- *Map witness test.* The _stronger_ condition $(I times.o Phi_W)(C) succ.eq.not 0$, evaluated as a smallest-eigenvalue check on the ampliation of the witness map; this can fire where the scalar test does not.
- *DPS relaxation.* The level-$2$ DPS relaxation; a feasible dual certificate flags entanglement.

The witness library and the DPS test are complementary: precomputed witnesses are cheap to evaluate but are individually much less likely to detect entanglement, whereas the DPS relaxation is more robust but far more expensive. Running all three and recording each score also lets us compare the methods directly, how often the witnesses alone suffice versus how much the DPS relaxation adds (section @results-section).

The components above assemble into a three-stage pipeline. The first two stages build two pools independently, a library of PNCP witness maps and a pool of bound entangled PPT states, and the third stage tests all of their compositions. Separating generation from testing lets each pool be built once, checkpointed, and reused across runs.

#algorithm-figure(
  "PPT2 counterexample search",
  {
    Comment[Stage 1: witness library]
    For($i = 1, ..., N_W$, {
      Line[Construct a PNCP map via @kmsz, rationalize its certificate (section @rationalization), and store the witness $W_i$.]
    })
    Comment[Stage 2: bound entangled candidate pool]
    For($j = 1, ..., N_C$, {
      Line[Sample a random PPT state via @ppt-gen; keep it only if we can verify it is entangled.]
    })
    Comment[Stage 3: test all compositions]
    For($"every ordered pair" (a, b)$, {
      Line[Form the composite channel $C = (I_n times.o Phi_a)(C_(Phi_b))$.]
      Line[Flag $C$ if any witness fires ($tr[W_i C] < 0$ or $(I times.o Phi_(W_i))(C) succ.eq.not 0$) or the level-$2$ DPS relaxation reports entanglement.]
      Line[If $C$ is entangled, re-validate exactly and export the certificate: a counterexample to the PPT2 conjecture.]
    })
  },
) <pipeline>

=== Direct PPT2 search by see-saw <gen-witness-ppt2>

The pipeline above tests a finite pool of compositions. We can instead search the composition manifold directly, asking whether a witness $W$ can be made to fire anywhere on it. Sharpening @min-ppt-witness, we minimise $W$'s expectation not over the whole PPT cone but over states that are themselves compositions of two PPT maps:
$ "minimize" quad & tr[W dot C_(rho_1 compose rho_2)] \ "subject to" quad & rho_(1,2) succ.eq 0, \ & rho_(1, 2)^Gamma succ.eq 0, \ & tr[rho_(1, 2)] = 1. $ <min-ppt2-witness>
A composition of PPT maps is itself PPT, so a negative optimum would exhibit a PPT _and_ entangled composition, a PPT2 counterexample witnessed by $W$. The ampliation is bilinear in the pair $(rho_1, rho_2)$, so @min-ppt2-witness is not a single SDP but a bilinear matrix inequality; we solve it by _see-saw_, freezing one factor and optimising the other in alternation, restarting from several random PPT initialisations. Being non-convex, see-saw reaches only a local optimum, so a non-negative result is supporting evidence rather than proof; the PPT-cone problem @min-ppt-witness is its convex relaxation and lower bound.

=== A note on the representation freedom <asym-witnesses>

As a by-product we explored the representation freedom of section @gram-freedom. Each witness can be expanded into a family of Gram representatives $M_0 + sum_alpha lambda_alpha N_alpha$ over a basis ${N_alpha}$ of $L$. We do not, in the end, rely on this. The relations spanning $L$ vanish on _real_ product vectors only, so an asymmetric representative stays block positive over $RR$ but need not remain so over $CC$. Guaranteeing such behavior is to our knowledge also an NP-hard problem without a simple solution in sight. The canonical, IPT representative, annihilating the skew-symmetric part as previously proposed remains the preferred approachfor our search, together with generic, non-symmetrized candidates. We leave a principled complex-domain extension to future work.

== Implementation Architecture

The implementation is a small Julia package in the `code/` directory. The core logic lives in the `ppt2` module (`code/src/ppt2.jl`), with the PNCP construction split into `code/src/pncp.jl` and included into the same module. Together they map one-to-one onto the operations defined above:

#table(
  columns: (1fr, 1.8fr),
  align: (left, left),
  stroke: none,
  table.header([*Function*], [*Role*]),
  table.hline(),
  [`rand_ppt`], [Sample a random PPT state (@ppt-gen): the positivity shift, with optional block symmetrization (`ppt_invariant`).],
  [`ampliation`], [Compute $(I times.o Phi)(C)$, used both to compose maps and to apply a witness test.],
  [`sample_pncp_form`, `segre_kernel_basis`, `non_sos_form`], [The KMSZ construction: sample Segre-variety points, build the linear forms $h_j$, and produce the non-SOS quadratic form $f$.],
  [`solve_sos`], [Set up and solve the SOS feasibility/optimization SDP for a given relaxation degree $l$; optionally trigger rationalization.],
  [`rationalize_certificate`], [Post-solver rationalization: zero the first $e$ Gram eigenvalues, recover rational coefficients, and re-check non-SOS.],
  [`find_pncp_poly`, `pncp_mat`, `poly2mat`], [Orchestrate witness generation with retries and export the certificate as a Choi matrix.],
  [`min_ppt_witness`, `min_ppt2_witness`], [The witness-restricted SDPs @min-ppt-witness and @min-ppt2-witness: bound entangled states from a witness, and the see-saw search over compositions.],
)

The module also exports a few supporting primitives: `rand_sep` and `rand_psd` for reference states, `is_ppt` for the PPT check, `gram_freedom` and `is_block_positive` for the witness-representation freedom and the block-positivity check of section @gram-freedom, and `swap`/`antisymmetric_projector` for the antisymmetric-subspace construction of @Sindici_2018.

The package leans on the established Julia optimization and quantum-information stack rather than reimplementing it. Polynomials and the SOS cone are handled by `DynamicPolynomials` and `SumOfSquares`; the resulting semidefinite programs are modelled with `JuMP` and solved by `MOSEK` through `MosekTools` (any JuMP-compatible SDP solver could be substituted). The DPS hierarchy is not reimplemented: the search driver calls `entanglement_robustness` from `Ket`, an existing quantum-information toolbox that also supplies utilities such as the partial transpose. Matrices and witness libraries are serialized with `JLD2` so that generation and search can be separated and resumed.

Several command-line drivers in `code/scripts/` orchestrate the long-running jobs, all sharing a `common.jl` harness that provides resumable, reproducible, multithreaded batch generation: completed batches are detected and skipped on a rerun, and every candidate is seeded deterministically so a configuration yields the same dataset regardless of thread count. The generation drivers build the two pools: `gen_pncp.jl` constructs the PNCP witness library (Stage 1), while `gen_ppt.jl` samples and DPS-filters a random bound entangled candidate pool and `gen_witness_ppt.jl` produces, for each witness, the bound entangled state the witness-restricted SDP extracts from the full PPT cone (Stage 2). `test_ppt2.jl` then runs the threaded all-pairs search of @pipeline (Stage 3), logging any detection together with the offending state and witness. The remaining drivers support the analysis of section @results-section: `compare_detection.jl` records every criterion's score on each detected state, the two witness tests and the DPS robustness, so the detection methods can be compared directly; `gen_witness_ppt2.jl` sharpens the witness-restricted SDP from the whole PPT cone down to the composition manifold by a see-saw over PPT-map pairs; and `cross_trace.jl` and `cross_ampl.jl` measure how broadly each witness detects states beyond its own. The `code/test/` suite checks the construction against reference values and verifies that generated maps are positive on large random samples, and the `code/notebooks/` directory documents the rationalization, PPT-state, and UPB workflows interactively.

= Results <results-section>

This chapter summarizes computational results from the pipeline implementations.

== Search scale and the broad composition search <broad-search>

All experiments are in the $4 times 4$ case, the smallest dimension where the PPT2 conjecture remains open; the methods extend to slightly larger dimensions, but the cost grows quickly and we did not pursue them.

Stage 1 produced a library of #nbr[10,000] PNCP witness maps via the KMSZ construction (@kmsz), each verified positive but not completely positive through the rationalization and SOS checks of section @rationalization. Stage 2 produced the bound entangled candidate pool in two ways: random sampling with a level-2 DPS entanglement filter (@ppt-gen) kept #nbr[849] of #nbr[1,435] sampled states (a hit rate of #nbr[59%]), and the witness-restricted SDP @min-ppt-witness extracted a further #nbr[10,000] bound entangled states, one per detecting witness. Building the two pools took approximately #nbr[5] hours on #nbr[80] CPU cores.

Stage 3 screens compositions by the witness criteria of @pipeline -- the scalar and map scans -- reserving the heavier level-2 DPS relaxation for any flagged pair. Testing a single pair takes on the order of seconds even when multithreaded, so the full quadratic sweep over a pool of thousands of states is infeasible; we therefore scan as large a slice as the budget allows. The portion on record covers #nbr[5,000] ordered pairs, with #nbr[zero] detections. This pairwise sweep is necessarily partial, which is exactly why the witness-restricted SDP probe of the next section -- complete over the entire witness library -- carries the weight of the argument.

No composition was found entangled in any test: every composite of two PPT maps screened as PPT and separable to within tolerance, leaving the PPT2 conjecture without a counterexample.

== The witness-restricted SDP probe <ppt2-witness-probe>

The broad search asks whether any composition in a finite pool is entangled. A sharper question is whether the witnesses we hold can detect entanglement on the composition manifold _at all_, even when the state is chosen to make them fire as hard as possible. The two witness-restricted SDPs of Sections @gen-witness-ppt and @gen-witness-ppt2 answer it.

The first, `min_ppt_witness` (@min-ppt-witness), minimises $tr[W rho]$ over the entire PPT cone. A negative optimum exhibits a PPT entangled state that $W$ detects. Over the library, #nbr[all 10,000] witnesses attained a negative optimum, with margins ranging from about #nbr[$10^(-7)$] to #nbr[$3 times 10^(-2)$] (median #nbr[$5 times 10^(-5)$]). This confirms that the canonical witnesses do carry genuine PPT-detection power once the state is optimised to them, as Theorem @indecomp-thm predicts, in contrast to their inertness on randomly drawn states (section @witness-findings); the small typical margin is itself the near-boundary effect.

The second, `min_ppt2_witness` (@min-ppt2-witness), restricts that minimisation from the whole PPT cone down to the composition manifold via the see-saw of section @gen-witness-ppt2, here run from #nbr[32] random restarts of at most #nbr[64] steps each.

The contrast between the two is the central result of this chapter. Every one of the 10,000 witnesses is live over the full PPT cone; yet when we ran the heavier see-saw on a sample of #nbr[2,000] of them, #nbr[not one] detected anything on the composition manifold. Every see-saw optimum was non-negative, lying in #nbr[$[2 times 10^(-16), 2 times 10^(-9)]$] — flush against the separable boundary rather than crossing it. The witnesses that are demonstrably _live_ over the PPT cone go _dead_ once the state is constrained to be a composition of two PPT maps, exactly the signature expected if PPT2 holds in $4 times 4$: the composition manifold lies within the entanglement-breaking region, beyond the reach of every witness we can construct.

== Detection methods compared <method-comparison>

Recording every criterion's score on each state lets us compare the two detection routes directly, the comparison made meaningful by the DPS--PNCP correspondence. In practice they are far from equivalent. Of the #nbr[849] random bound entangled states that the level-2 DPS relaxation certifies (section @broad-search), the PNCP witness library detects #nbr[none] by the scalar test: all #nbr[8.5 million] witness--state evaluations were strictly positive, the smallest being #nbr[$+148$]. The tight, near-decomposable KMSZ witnesses thus entirely miss the generic bound entangled states that the DPS relaxation reaches, the inertness of section @witness-findings made quantitative.

Within the witness scan itself, the map witness test is the stronger of its two forms: the inequality $lambda_min ((I times.o Phi^*)(rho)) <= tr[C_Phi rho]$ of section @qse-section guarantees it detects everything the scalar test does, and in general strictly more, since it inspects the whole spectrum of $(I times.o Phi_W)(C)$ rather than a single expectation.

== Witness generality: cross-detection <cross-detection>

Each witness is built tight to one prescribed family of zeros, so it is natural to ask how far it reaches beyond the single state it was made for. We cross-evaluated every witness against the bound entangled state that the witness-restricted SDP extracted for every _other_ witness, separately for the scalar criterion (`cross_trace.jl`) and the map criterion (`cross_ampl.jl`). The result is stark: across all #nbr[$10^8$] witness--state pairs, neither criterion produced a single foreign detection. Both tests flagged only the #nbr[10,000] diagonal pairs, each witness detecting exactly the one state it generated and #nbr[no] other. This vanishing cross-detection rate is the sharpest expression of the tightness of the witnesses: each cuts the entangled region in just one place, at its own prescribed state, so the library behaves as #nbr[10,000] essentially disjoint detectors rather than a redundant cover.

== Representation and symmetry effects <witness-findings>

Two further observations concern the witness representations themselves rather than the search outcome.

First, the partial-transpose-invariant (symmetric) representatives are inert for our purpose. Symmetrized witnesses, the $M_0 = M_0^Gamma$ representatives returned by the construction, never flagged a single PPT entangled state in our experiments, and PPT states built in the symmetric `ppt_invariant` shape never appeared as counterexamples. This is _not_ because they are decomposable: the canonical KMSZ witness is exactly this partial-transpose-invariant representative (section @complexification), and it is provably indecomposable (Theorem @indecomp-thm), so it does carry genuine PPT-detection power. The inertness is instead a near-boundary effect. As argued previously, these witnesses are only tight and sit close to the decomposable cone, so whatever PPT entanglement they detect does so with a margin at the level of solver noise -- consistent with the small PPT-violation margins (of order $10^(-6)$) reported for the same construction by @Masse_2026. In practice such detections fall below our tolerance, so a fixed canonical witness rarely fires on a randomly drawn bound entangled state. The robust detection in the all-pairs search therefore rests on the DPS relaxation; the witness library's genuine reach is exposed instead when the state is optimised to the witness, which is the point of the witness-restricted SDP probe (section @ppt2-witness-probe). For the cheap per-state scan we use the canonical representative on generic, non-symmetrized candidates.

Second, we also probed the asymmetric representatives directly (section @asym-witnesses). On the real domain they behave exactly as the theory predicts: they remain positive on all real product vectors, and over many random trials never assigned a negative score to a real separable state. Whether this extends to the complex domain is open: the relations underlying the representation freedom (section @gram-freedom) vanish only on real product vectors, so an asymmetric representative that is a valid witness over $RR$ need not be block-positive over $CC$. Without an analytic guarantee the asymmetric witnesses proved too unreliable to use: the map test in particular can return a negative eigenvalue on a separable input -- a spurious detection -- precisely because an asymmetric representative need not be a positive map over $CC$. Establishing complex block-positivity, or restoring it by a suitable projection, is left for future work (section @future-work).

A more direct comparison may also be with the recent paper @Masse_2026, where the authors produce a library of 20000 KMSZ witnesses in $3 times 3$, stating it took "several days", albeit on a marginally weaker machine. Our result is still orders of magnitude faster.

== Performance

A direct head-to-head comparison with the earlier MATLAB prototype is difficult, because the two implementations were built with different goals and make different design choices. The prototype drew random _integer_ matrices to produce cosmetically nicer certificates; this is not random enough and frequently failed, so a single candidate often needed many recomputation attempts. Our implementation can still generate candidates that way, but for speed it samples from a normal distribution, which almost always succeeds on the first try. The two also rationalize differently: the prototype rationalizes the entire problem up front, paying for rational arithmetic throughout, whereas we rationalize only after the SDP is solved (@witness-library). The solver itself is largely language-independent, but the Julia code benefits from more efficient data handling and tighter solver integration.

With those caveats, we report rough wall-clock figures for the current implementation. Generating both pools for the $4 times 4$ search -- the library of 10,000 PNCP witness maps and the bound entangled candidate pool -- took approximately #nbr[5] hours on #nbr[80] CPU cores, i.e. on the order of #nbr[?? s per PNCP map -- need per-map timing] per map amortised. For the same dimensions the prototype took orders of magnitude longer to produce candidates. #margin-note[fill the MATLAB-prototype figure (measured or from the article) and the per-map Julia timing]

= Discussion <discussion>

This thesis develops a computational approach to the PPT2 conjecture and applies it to the smallest open case, $4 times 4$. The pipeline generates PPT candidate states, builds a library of provably indecomposable PNCP witnesses from the KMSZ construction, and screens composite channels for entanglement by two distinct routes -- entanglement witnesses, and the level-2 DPS hierarchy -- with a post-solver rationalization step that turns any numerical detection into an exact certificate.

No composition was ever found to be entangled. The strongest evidence comes not from the brute pairwise sweep, which the per-pair cost keeps necessarily partial, but from the witness-restricted SDP probe (section @ppt2-witness-probe): although every one of our witnesses detects a bound entangled state somewhere in the full PPT cone, not one detects anything once the state is constrained to be a composition of two PPT maps, with the see-saw optima sitting exactly on the separable boundary. Witnesses that are demonstrably live over the PPT cone go uniformly dead on the composition manifold, precisely the signature expected if the PPT2 conjecture holds.

This is computational evidence, not a proof, and it inherits the limitations of its tools. The KMSZ witnesses, though provably indecomposable (Theorem @indecomp-thm), are only tight and lie close to the decomposable cone, so they detect PPT entanglement weakly. Much of the detection burden fell on the DPS relaxation, itself confined to level 2 in this dimension. A genuine counterexample, should one exist, might also lie in a region that our random candidate generation never reaches. What the methodology does provide is a reproducible, exactly-validated framework that returns a clean negative across every test we could run, and that extends in principle to the next open dimensions.

== Complexity and Practical Limits <complexity-limits>

All the problems we are looking at are SDPs with exponential @Gharibian_2009 complexity in dimension and relaxation depth @Doherty_2004. Practical bottlenecks include:
1. processing power: these problems are by nature not easily parallelizable, simplex methods are intrinsically sequential and internal point methods rely on factoring large _sparse_ matrices, so the potential of GPU speedups is fairly low @MOSEK.
2. memory growth: the size of the SDP grows exponentially with dimension and relaxation depth, leading to memory bottlenecks even for moderate dimensions, i.e. for $4 times 4$ states DPS level 3 is already infeasible on standard hardware, requiring hundreds of gigabytes of RAM. Even improvements such as adding KKT constraints to achieve finite convergence @Harrow_2017 lead to significant increases in problem size, so they are not without cost.
3. solver instability near degeneracy: SDP solvers produce floating-point solutions, and when the problem is near the boundary of feasibility, small numerical errors can lead to incorrect conclusions about positivity. This is particularly problematic for our purposes, since we are interested in the very existence of such bound states.
4. The generated witnesses are provably indecomposable (Theorem @indecomp-thm), so each one can in principle detect some PPT entangled state. They are, however, only tight and therefore lie close to the decomposable cone, detecting PPT entanglement only weakly. This compounds the numerical fragility of point 3, and how best to improve a given witness, remains open.
5. Searching for PPT candidates in the first place is non-trivial, and random generation may not be sufficient to find counterexamples if they exist. We may require a more structured approach if the volume of the search space is in fact 0 (this would not mean the conjecture holds).

== Future work <future-work>

Several directions remain open.

The most immediate is the _complex-domain extension of the representation freedom_ (Sections @gram-freedom, @asym-witnesses). The asymmetric Gram representatives are valid witnesses over $RR$ but lose block-positivity over $CC$, because the Segre relations spanning $L$ vanish only on real product vectors. Klep et al. restore positivity over $CC$ for the canonical representative by annihilating the skew-symmetric part @Klep_2017; an analogous analytic guarantee for the asymmetric family, or a projection that restores complex block-positivity while keeping its sharper cut of the entangled region, would turn the freedom from a real-only curiosity into a usable supply of independent witnesses. We pursued this only numerically and found the evidence too weak to rely on, so the analytic question is left open.

A second direction is _structured candidate generation_. Random PPT sampling discards most of its draws as separable and leaves unstructured survivors (limitation 5, section @complexity-limits); if the volume of bound entangled compositions is in fact vanishing, a random search will never reach a counterexample even if one exists. The constructions surveyed in section @ppt-gen -- unextendible product bases @Bennett_1999, antisymmetric-subspace states @Sindici_2018, and symmetric random induced states @Louvet_2025 -- offer entanglement-by-design pools that a targeted search could exploit.

Finally, the pipeline is dimension-limited by the cost of the SDPs (section @complexity-limits). Pushing beyond $4 times 4$, whether by exploiting state symmetry @Britz_2025, by the finite-convergence KKT extensions of @Harrow_2017, or by specialised solvers, would let the same methodology test the conjecture where even less is known.

#bibliography(
  title: "References",
  style: "ieee",
  "bibliography.bib"
)
