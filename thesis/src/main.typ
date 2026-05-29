#import "conf.typ": conf
#import "@preview/drafting:0.2.2": *
#import "@preview/algorithmic:1.0.7"
#import "@preview/cheq:0.3.0": checklist
#import algorithmic: *
#import "@preview/lemmify:0.1.8": *

#show: style-algorithm
#show: checklist

#show: conf.with(
  title_en: "A Software Approach to the PPT2 Conjecture",
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

#let (theorem, definition, rules: thm-rules) = default-theorems("thm-group", lang: "en", thm-numbering: thm-numbering-linear)
#show: thm-rules

// fixes "layout did not converge" warning
// #let margin-note = margin-note.with(dy: 0pt)
#let reference = margin-note("reference")
#let PPT2 = $"PPT"^2$

= TODO <todo>
*Writing*
- [ ] Write the abstracts: English, Slovene, and the extended Slovene summary.
- [ ] Polish the Introduction: spell out the implications of the conjecture holding versus failing, and sharpen the statement of current knowledge ($n <= 3$ proven, $n >= 4$ open).
- [ ] Expand the Conclusion into a genuine discussion: the Gram-representation freedom, the inertness of the symmetric (partial transpose invariant) family, the limitations, and the outlook.

*Results* (pending)
- [ ] Finalize the PPT2 search results (Section @results-section): final candidates and witnesses, hardware, wall-clock, and outcome.
- [ ] Fill in the performance numbers in the Performance section: Julia generation timings and the MATLAB-prototype comparison.

*Figures and references*
- [ ] Add the flagged figures: the separable set with a witness hyperplane, maybe a schematic of the pipeline.
- [ ] Complete and verify all citations: resolve the remaining reference placeholders and margin notes.

*Open questions / future work*
- [ ] Investigate the complex-domain extension of the asymmetric witnesses (Sections @gram-freedom, @witness-findings): block-positivity over $CC$ fails, how to fix it.
- [ ] Consider structured candidate generation (UPB or symmetric random induced states) should the random search volume prove effectively zero (cf. limitation 5 in @complexity-limits).

= Introduction <intro>
Quantum entanglement is a central nonclassical resource in quantum information theory. It supports communication and cryptographic tasks and provides a structural lens for understanding quantum channels @Horodecki_2009.

This thesis studies the #PPT2 conjecture: whether the composition of two PPT maps is always entanglement breaking @Christandl_2019. The conjecture has a clear operational interpretation, since entanglement-breaking maps destroy all bipartite entanglement with any reference system.

Known progress motivates a computation-first strategy. The conjecture is proven in low-dimensional regimes and in specific families such as Choi-type maps @Chen_2019 @Singh_2022. In higher dimensions, the geometry of positive maps becomes more complex, and direct analytic classification remains difficult.

The core challenge is twofold. First, separability testing is hard in the regimes where the conjecture remains open. Second, numerical certificates near feasibility boundaries are fragile. Therefore, this thesis focuses on methodology: how to construct, compose, test, and validate candidate objects in a reproducible way.

The main objectives are:
1. Design a Julia workflow for candidate generation, map composition, and entanglement checks.
2. Integrate DPS-based semidefinite relaxations as a baseline witness route.
3. Integrate an SOS-based PNCP witness-construction route with post-solver validation.

The rest of the thesis is organized as follows. Chapter 2 provides mathematical prerequisites and formulations. Chapter 3 presents the computational workflow and implementation architecture. Chapter 4 discusses methodological limitations and numerical reliability. Chapter 5 concludes with implications and next directions.

= Theoretical Background
This chapter establishes the mathematical framework underpinning the PPT2 conjecture and its computational study. We fix notation and define the key objects: linear maps, their positivity properties, quantum states, and entanglement criteria, then develop the semidefinite programming tools used in the implementation.

== Preliminaries

=== Notation and matrix spaces

We write $M_n$ for the algebra of complex $n times n$ matrices, $H_n subset.eq M_n$ for the subspace of Hermitian matrices ($A = A^*$), and $M_n^+ subset.eq H_n$ for the cone of positive semidefinite matrices.

#definition(name: "Positive semidefinite and positive definite")[
  A matrix $A in H_n$ is _positive semidefinite_ (PSD), written $A succ.eq 0$, if $x^* A x >= 0$ for all $x in CC^n$, or equivalently if all eigenvalues of $A$ are non-negative. If all eigenvalues are strictly positive, $A$ is _positive definite_ (PD), written $A succ 0$.
]

The _Hilbert-Schmidt inner product_ on $H_n$ is
$ chevron.l A, B chevron.r := "Tr"[A B], $
making $H_n$ an inner product space. For general $M_n$ the inner product is $chevron.l A, B chevron.r = "Tr"[A^* B]$.

We work with bipartite systems on the tensor product space $CC^m times.o CC^n$ (subsystem dimensions $m$ and $n$). The operator space is $M_m times.o M_n tilde.eq M_(m n)$. We write $E_(i j) in M_n$ for the matrix unit with a $1$ in position $(i,j)$ and zeros elsewhere.

=== Linear maps on matrix spaces

#definition(name: "Positive, k-positive, and completely positive maps")[
  Let $Phi: M_n -> M_m$ be a linear map.
  - $Phi$ is _positive_ (P) if $A succ.eq 0 => Phi(A) succ.eq 0$.
  - $Phi$ is _$k$-positive_ if $I_k times.o Phi: M_k times.o M_n -> M_k times.o M_m$ is positive.
  - $Phi$ is _completely positive_ (CP) if it is $k$-positive for every $k in NN$ @Chen_2019.
]

Every CP map is positive. Maps that are positive but not completely positive (PNCP) play a central role: they detect entanglement in a way that CP maps cannot.

#definition(name: "Choi-Jamiolkowski isomorphism")[
  The _Choi matrix_ of a linear map $Phi: M_n -> M_m$ is
  $ C_Phi = sum_(i,j=1)^n E_(i j) times.o Phi(E_(i j)) in M_n times.o M_m. $
  The assignment $Phi arrow.bar C_Phi$ is a linear isomorphism between maps $Phi: M_n -> M_m$ and matrices in $M_n times.o M_m$. The map is recovered by $Phi(A) = "Tr"_1 [(A^T times.o I_m) C_Phi]$.
]

Under this isomorphism, $Phi$ is CP if and only if $C_Phi succ.eq 0$ @Choi_1975, and $Phi$ is trace-preserving if and only if $"Tr"_B [C_Phi] = I_n$.

#definition(name: "Partial transpose")[
  For a product matrix $A times.o B in M_m times.o M_n$, the _partial transpose_ with respect to subsystem $B$ is defined as
  $ (A times.o B)^(Gamma_B) := A times.o B^T. $
  Extending by linearity to a general $rho in M_m times.o M_n$, written in block form as $rho = (rho_(i j))_(i,j=1)^m$ with blocks $rho_(i j) in M_n$:
  $ rho^Gamma := (I_m times.o T)(rho) = (rho_(i j)^T)_(i,j=1)^m, $
  where $T$ is the transposition map. For Hermitian $rho$, the partial transpose with respect to either subsystem yields the same spectrum: $(rho^(Gamma_A))^T = rho^(Gamma_B)$, so $rho^(Gamma_A) succ.eq 0 <=> rho^(Gamma_B) succ.eq 0$. We therefore write $rho^Gamma$ without specifying the subsystem.
]

#definition(name: "PPT state and PPT map")[
  A state $rho in M_m times.o M_n$ is _PPT_ if $rho succ.eq 0$ and $rho^Gamma succ.eq 0$. A map $Phi: M_n -> M_m$ is _PPT_ if its Choi matrix is PPT:
  $ C_Phi succ.eq 0 quad "and" quad C_Phi^Gamma succ.eq 0. $
]

Since $C_Phi succ.eq 0$ characterizes CP maps, a PPT map is a CP map whose Choi matrix has non-negative partial transpose. The composition of two PPT maps is trivially PPT.

=== Map composition via Choi matrices <map-comp>

The Choi matrix of the composition $Phi compose Psi$ is not the product $C_Phi C_Psi$. Instead, it is obtained by applying the _ampliation_ $I_n times.o Phi$ to $C_Psi$:
$ C_(Phi compose Psi) = (I_n times.o Phi)(C_Psi), $
where the right-hand side applies $Phi$ block-wise to each $n times n$ block of $C_Psi$. Expanding in the standard basis yields the index formula
$ (C_(Phi compose Psi))_(i p, j q) = sum_(k,l) (C_Psi)_(i k, j l) (C_Phi)_(k p, l q). $
This formula is used directly in the implementation to compose two candidate PPT maps from their Choi matrices.

== PPT2 Conjecture

With PPT maps and entanglement-breaking maps (Definition @eb-def below) in hand, we can now state the central object of study.

#theorem(name: "PPT2 Conjecture")[
  If $Phi_1$ and $Phi_2$ are PPT maps, then $Phi_1 compose Phi_2$ is entanglement breaking @Christandl_2019.
]<ppt2>

The conjecture is proven for $n = 2$ (since all PPT states in $M_2 times.o M_n$ are separable @Horodecki_1996), for $n = 3$ @Chen_2019, and for Choi-type maps in all dimensions @Singh_2022. It remains open for $n >= 4$ in the general case. The rest of this chapter develops the formalism needed to study the conjecture computationally.

== Quantum States and Entanglement <qse-section>

A _quantum state_ on $CC^n$ is a density matrix: a PSD matrix $rho in M_n^+$ with $"Tr"[rho] = 1$. A _quantum channel_ is a completely positive trace-preserving (CPTP) map $Phi: M_n -> M_m$.

#definition(name: "Separability and entanglement")[
  A bipartite state $rho in M_m times.o M_n$ is _separable_ if
  $ rho = sum_i p_i rho_i^A times.o rho_i^B, quad p_i >= 0, quad sum_i p_i = 1, quad rho_i^A in M_m^+, quad rho_i^B in M_n^+. $
  Otherwise $rho$ is _entangled_.
]

Deciding separability is NP-hard in general @Gharibian_2009. Since an exact test is computationally intractable, one relies on necessary conditions: a state that fails such a condition must be entangled.

#definition(name: "Entanglement-breaking map")[
  A map $Phi: M_n -> M_m$ is _entanglement breaking_ (EB) if $(I_k times.o Phi)(rho)$ is separable for every $k in NN$ and every state $rho in M_k times.o M_n$ @Horodecki_2009.
] <eb-def>

Testing $k = 1$ alone would only check that $Phi$ preserves separability of bipartite states; the ampliation over all $k$ ensures $Phi$ destroys entanglement with any external reference system. Equivalently, $Phi$ is EB if and only if its Choi matrix $C_Phi$ is a separable state @Horodecki_2009.

=== Separability criteria

Since testing separability exactly is NP-hard, one uses necessary conditions. A state that violates such a condition is certifiably entangled; satisfying all known criteria does not guarantee separability. For a comprehensive survey see @Guhne_2009.

*PPT criterion.* If $rho$ is separable, then $rho^Gamma succ.eq 0$ @Peres_1996 @Horodecki_1996. More precisely: the transposition map $T$ is positive but not completely positive, so applying $I times.o T$ to an entangled state can produce a non-PSD result. If $rho^Gamma succ.eq.not 0$, then $rho$ is entangled. The condition is necessary but not sufficient: bound entangled (PPT entangled) states exist whenever $m n > 6$ @Horodecki_2009.

*Range criterion.* If $rho in M_m times.o M_n$ is separable, there exist product vectors $lr({|psi_i chevron.r times.o |phi_i chevron.r})$ spanning the range of $rho$ such that $lr({|psi_i chevron.r times.o overline(|phi_i chevron.r)})$ spans the range of $rho^Gamma$ @Horodecki_1997. The range criterion can detect certain PPT entangled states, but fails when $rho$ is full rank (e.g., under noise), since then any set of vectors spans its range trivially.

*Reduction criterion.* If $rho in M_m times.o M_n$ is separable, then
$ rho_A times.o I_n - rho succ.eq 0 quad "and" quad I_m times.o rho_B - rho succ.eq 0, $
where $rho_A = "Tr"_B [rho]$ and $rho_B = "Tr"_A [rho]$ are the _reduced states_ (obtained by tracing out one subsystem) @Horodecki_1999.

*CCNR criterion.* Define the _realignment_ of $rho$ as the matrix $R(rho)$ with entries $(R(rho))_(i m+k,, j n+l) = rho_(i n+j,, k n+l)$. If $rho$ is separable, then $||R(rho)||_1 <= 1$, where $||dot||_1$ is the trace norm @Chen_2003. The CCNR criterion is independent of the PPT criterion and can detect some PPT entangled states that the other criteria miss.

*Majorization criterion.* If $rho$ is separable, then $lambda(rho) prec.eq lambda(rho_A)$ and $lambda(rho) prec.eq lambda(rho_B)$, where $lambda(dot)$ denotes the non-increasingly ordered eigenvalue vector @Nielsen_2001. Here $a prec.eq b$ denotes _majorization_: $sum_(i=1)^k a_i <= sum_(i=1)^k b_i$ for all $k$, with equality at $k = n$. (Note: $prec.eq$ here is majorization order, distinct from the positive-definiteness notation $succ.eq$.) This criterion follows from the reduction criterion @Nielsen_2001 and shares its limitations.

Beyond criterion-based tests, _algorithmic_ approaches reformulate separability as a convex optimization problem. The DPS hierarchy (Section @dps-section) is the most powerful general-purpose method. Criteria based on covariance matrices and entanglement witnesses (including Bell inequalities) cover complementary cases but are less directly applicable in our setting @Guhne_2009.

=== Entanglement witnesses

#definition(name: "Entanglement witness")[
  A Hermitian operator $W$ is an _entanglement witness_ if $"Tr"[W sigma] >= 0$ for all separable $sigma$, and $"Tr"[W rho] < 0$ for some entangled $rho$.
]

The set of separable states is convex and closed. By the Hahn-Banach separation theorem; any point outside a closed convex set is separated from it by a hyperplane, for every entangled $rho in.not "SEP"$ there exists a Hermitian operator $W$ (a hyperplane in $H_(m n)$) with $"Tr"[W rho] < 0$ and $"Tr"[W sigma] >= 0$ for all $sigma in "SEP"$ @Horodecki_2009. Geometrically, each witness $W$ cuts off a half-space that contains $rho$ but not $"SEP"$; the intersection of all such half-spaces recovers $"SEP"$ exactly. #margin-note[figure: separable set and witness hyperplane]

Witnesses are one-sided: $"Tr"[W rho] < 0$ certifies entanglement, but $"Tr"[W rho] >= 0$ does not certify separability.

#definition(name: "Decomposable witness")[
  A witness $W$ is _decomposable_ if $W = P + Q^Gamma$ for some $P, Q succ.eq 0$; otherwise it is _non-decomposable_ @Lewenstein_2000.
]

Decomposable witnesses cannot detect PPT entangled states: for any PPT $sigma$, $"Tr"[(P + Q^Gamma) sigma] = "Tr"[P sigma] + "Tr"[Q sigma^Gamma] >= 0$. Only non-decomposable witnesses are useful for our purposes.

Under the Choi-Jamiolkowski isomorphism, every entanglement witness $W = C_Phi$ corresponds to a PNCP map $Phi$, and vice versa. This gives two ways to test a state, and the map-based one is strictly stronger. The scalar test asks only whether $"Tr"[C_Phi rho] < 0$, a single separating hyperplane. The map test asks whether the ampliation $(I_k times.o Phi)(rho)$ has a negative eigenvalue, $(I_k times.o Phi)(rho) succ.eq.not 0$. The latter holds exactly when some vector $|psi chevron.r$ satisfies
$ chevron.l psi | (I times.o Phi)(rho) | psi chevron.r = "Tr"[(I times.o Phi^*)(|psi chevron.r chevron.l psi|) rho] < 0, $
where $Phi^*$ is the adjoint map (also positive). The map therefore encodes an entire family of scalar witnesses $W_psi = (I times.o Phi^*)(|psi chevron.r chevron.l psi|)$, one per output vector $|psi chevron.r$; the Choi-matrix test $"Tr"[C_Phi rho]$ is the single member obtained from the maximally entangled $|psi chevron.r$, whereas the eigenvalue check optimizes over all $|psi chevron.r$ at once. Hence
$ lambda_min ((I times.o Phi^*)(rho)) <= "Tr"[C_Phi rho], $
so the map detects every state the scalar test does, and in general strictly more @Horodecki_2009. The search uses this directly, screening every composite with both tests (Section @methods).

== Semidefinite Programming

A _semidefinite program_ (SDP) is a convex optimization problem in which a linear objective is minimized subject to a linear matrix inequality. The standard primal form is

$ "minimize" &quad c^T bold(x) \
  "subject to" &quad F(bold(x)) := F_0 + sum_(i=1)^n x_i F_i succ.eq 0, $ <sdp-primal>

where $bold(x) in RR^n$ is the optimization variable, $c in RR^n$ is the cost vector, and $F_0, F_1, ..., F_n in H_d$ are fixed Hermitian matrices. The associated _dual SDP_ is

$ "maximize" &quad -"Tr"[F_0 Z] \
  "subject to" &quad Z succ.eq 0, \
               &quad "Tr"[F_i Z] = c_i, quad i = 1, ..., n. $ <sdp-dual>

The _Slater condition_ for @sdp-primal requires the existence of a strictly feasible point: some $bold(x)$ with $F(bold(x)) succ 0$. When it holds for both primal and dual, _strong duality_ holds: primal and dual optima coincide. When $c = 0$, @sdp-primal is a _feasibility problem_. If infeasible, a dual feasible $Z succ.eq 0$ with $"Tr"[F_i Z] = 0$ and $"Tr"[F_0 Z] > 0$ certifies infeasibility.

=== Interior point methods

SDPs are solved in practice by _interior point methods_ (IPMs), which follow the _central path_; a smooth trajectory through the strict interior of the feasible region parameterized by a barrier coefficient $mu -> 0$. A standard IPM applied to a $d times d$ SDP with $n$ variables reaches $epsilon$-accuracy in $O(sqrt(d) log(1 slash epsilon))$ iterations, each requiring $O(n^2 d^2 + n d^3)$ operations @Vandenberghe_1996. This thesis uses MOSEK @MOSEK, a state-of-the-art IPM solver for semidefinite programs.

A property of IPMs important for our use: solutions lie in the _strict interior_ of the feasible region ($F(bold(x)) succ 0$), which provides numerical slack around the boundary. This is exploited in the post-solver rationalization (Section @rationalization). However, all solutions are floating-point approximations, and any conclusion near a feasibility boundary requires explicit validation.

=== Entanglement testing as an SDP

Separability testing is NP-hard @Gharibian_2009. One standard formulation maximizes the expectation of an observable $M$ over separable states:
$ h_("SEP")(M) = max {"Tr"(M rho) : rho in "SEP"}. $
The DPS hierarchy (Section @dps-section) relaxes this problem by replacing SEP with a sequence of increasingly tight SDP-representable outer approximations.

== DPS Hierarchy <dps-section>

The Doherty-Parrilo-Spedalieri (DPS) hierarchy @Doherty_2004 provides a sequence of SDP relaxations of separability that is _complete_: for every entangled state there exists a finite level $k$ at which the test detects it. In principle, running the hierarchy to convergence solves separability exactly; in practice, only the first few levels are computationally feasible.

=== Symmetric extensions

For a separable state $rho = sum_i lambda_i x_i x_i^* times.o y_i y_i^*$, define the $k$-fold _symmetric extension_
$ rho_k = sum_i lambda_i x_i x_i^* times.o (y_i y_i^*)^(times.o k) in H(CC^m times.o (CC^n)^(times.o k)). $
This extension satisfies three properties, each with a physical interpretation:

$ (I_m times.o Pi_k) rho_k (I_m times.o Pi_k) = rho_k, $ <ext-sym>
$ rho_k^(Gamma_s) succ.eq 0 quad forall s = 1,...,k, $ <ext-ppt>
$ "Tr"_(B_2 ... B_k)[rho_k] = rho. $ <ext-marg>

Here $Pi_k$ projects onto the bosonic (fully symmetric) subspace of $(CC^n)^(times.o k)$. Condition @ext-sym states that $rho_k$ is invariant under permutations of the $k$ copies of $B$. Condition @ext-ppt states that all partial transposes of $rho_k$ are PSD. Condition @ext-marg states that tracing out the extra $k-1$ copies of $B$ recovers $rho$.

#definition(name: "DPS set at level k")[
  $"DPS"_n^k$ is the set of states $rho in H(CC^m times.o CC^n)$ for which there exists a symmetric extension $rho_k in H(CC^m times.o (CC^n)^(times.o k))$ satisfying @ext-sym, @ext-ppt, and @ext-marg.
]

Each $"DPS"_n^k$ is defined by semidefinite constraints, so membership is testable in polynomial time for fixed $k$. The hierarchy satisfies @Doherty_2004:
1. $"SEP"_n subset.eq "DPS"_n^k$ and $"DPS"_n^(k+1) subset.eq "DPS"_n^k$ for all $k >= 1$.
2. $"DPS"_n^1$ is equivalent to the PPT criterion.
3. Asymptotic completeness: $inter.big_(k >= 1) "DPS"_n^k = "SEP"_n$.

=== Feasibility SDP and witness extraction

Testing $rho in "DPS"_n^k$ amounts to searching for $rho_k$ satisfying @ext-sym, @ext-ppt, and @ext-marg: this is an SDP feasibility problem. When it is infeasible (no valid extension exists), the dual variable $Z succ.eq 0$ yields an entanglement witness.

Concretely, the SDP dual to the level-$k$ test has a feasible $Z$ whenever $rho$ is entangled at that level. The operator
$ W = "Tr"_(B_2 ... B_k) [Z] $
is an entanglement witness for $rho$ @Doherty_2004: one can verify $"Tr"[W rho] < 0$ directly from $Z$, without relying on the floating-point primal solution. This dual certificate is how we confirm any entanglement detection numerically.

=== Improvements and limitations

Several enhancements to the basic hierarchy are known. Harrow, Natarajan, and Wu @Harrow_2017 add first-order optimality (KKT) conditions to the DPS SDP, achieving _finite convergence_: infeasibility is certified at a finite level rather than only asymptotically. The KKT conditions are linear constraints on the Lagrange multipliers of the original SDP; they increase the variable count substantially and make clean witness extraction from the dual more involved.

Specialized hierarchies have been developed for states with symmetry. For diagonal unitary invariant states, Britz and Laurent @Britz_2025 give a drastically smaller SDP at each level. For Werner states and isotropic states, explicit separability conditions are known. However, none of these apply to our search: the PPT2 conjecture is already proven for the relevant symmetric families @Singh_2022, so our search must use the general DPS hierarchy.

== Sum-of-Squares Polynomials and PNCP Maps

A complementary approach to DPS exploits the polynomial representation of linear maps. The key insight: positivity of a map is equivalent to non-negativity of a polynomial, and complete positivity is equivalent to that polynomial being a sum of squares (SOS). Since SOS is characterized by a semidefinite condition (via Gram matrices), this connects map positivity directly to SDPs.

=== Non-negative polynomials and SOS

Let $RR[bold(x), bold(y)]$ be the ring of real polynomials in $bold(x) in RR^n$, $bold(y) in RR^m$. A polynomial $p$ is _non-negative_ if $p(bold(x), bold(y)) >= 0$ for all real inputs; it is a _sum of squares_ (SOS) if $p = sum_i q_i^2$ for polynomials $q_i$. Every SOS polynomial is non-negative, but not every non-negative polynomial is SOS. The gap between the two is the source of PNCP maps.

=== Gram matrix and SDP connection

#definition(name: "Gram matrix representation")[
  A homogeneous polynomial $p$ of degree $2d$ is SOS if and only if there exists $G succ.eq 0$ such that
  $ p(bold(z)) = bold(v)(bold(z))^T G bold(v)(bold(z)), $ <gram-rep>
  where $bold(v)(bold(z))$ is the vector of monomials of degree $d$ in $bold(z)$.
]

Testing whether $p$ is SOS therefore reduces to a semidefinite feasibility problem: find $G succ.eq 0$ satisfying the linear constraints that equate coefficients of @gram-rep with those of $p$. If no such $G$ exists, $p$ is not SOS.

=== Polynomial representation of maps

Each linear map $Phi: M_n -> M_m$ corresponds to a biquadratic polynomial:
$ p_Phi (bold(x), bold(y)) := bold(y)^T Phi(bold(x) bold(x)^T) bold(y). $

"Non-negative" here means $p_Phi (bold(x), bold(y)) >= 0$ for all real $(bold(x), bold(y))$, the polynomial analogue of positivity for matrices. The fundamental correspondence @Klep_2017 is:
- $Phi$ is _positive_ if and only if $p_Phi$ is non-negative on $RR^n times RR^m$.
- $Phi$ is _completely positive_ if and only if $p_Phi$ is SOS.

A PNCP map corresponds exactly to a non-negative non-SOS polynomial, and constructing such a polynomial gives an entanglement witness directly.

=== Non-uniqueness of the matrix representation <gram-freedom>

The polynomial $p_Phi$ does not determine the matrix $C_Phi$ uniquely. Writing $p_Phi$ as a quadratic form in the product monomials $bold(z) = bold(x) times.o bold(y)$,
$ p_Phi (bold(x), bold(y)) = bold(z)^T M bold(z), quad bold(z) = bold(x) times.o bold(y), $
any symmetric $M$ that reproduces the coefficients of $p_Phi$ is an admissible representation. Two such matrices $M$, $M'$ give the same polynomial precisely when $M - M'$ vanishes on all product vectors, $bold(z)^T (M - M') bold(z) = 0$. The matrices with this property form a linear space $L$, spanned by the $2 times 2$ minor (Segre) relations $(x_i y_k)(x_j y_l) - (x_i y_l)(x_j y_k)$; concretely $L tilde.eq and^2 RR^n times.o and^2 RR^m$, with $dim L = binom(n, 2) binom(m, 2)$. The admissible representations therefore form an affine space $M_0 + L$.

Every representative encodes the same polynomial, hence the same values on product vectors and, by linearity, the same expectation $"Tr"[M sigma]$ on every separable $sigma$. They differ only off the Segre variety, i.e. on entangled states. A single non-negative non-SOS polynomial thus yields not one witness but a whole family of them, all sharing the same separable boundary while cutting the entangled region differently. One subtlety carries into detection: the relations defining $L$ vanish on _real_ product vectors but not on complex ones, so members of $M_0 + L$ agree as witnesses over $RR$ but need not agree, or even remain block-positive, over $CC$.

#algorithm-figure(
  "KMSZ construction for PNCP maps",
  {
    Line[Sample random rational points $x^((1)), ..., x^((t)) in QQ^n$ and $y^((1)), ..., y^((t)) in QQ^m$.]
    Line[Form rank-1 bilinear forms $h_j (bold(x), bold(y)) = chevron.l x^((j)), bold(x) chevron.r dot chevron.l y^((j)), bold(y) chevron.r$. Each $h_j$ is a product of two linear forms, so $sum_j h_j^2$ is SOS.]
    Line[Find $f in.not "span"{h_1, ..., h_t}$ over $QQ$ such that $f$ is non-negative and not SOS. This uses linear algebra: choose $f$ orthogonal to the SOS cone within the space of degree-$2$ biquadratic forms @Klep_2017.]
    Line[Solve the SDP: find the smallest $delta > 0$ such that $F_delta := delta f + sum_j h_j^2 >= 0$. For small enough $delta$, the SOS term $sum h_j^2$ dominates $delta f$, ensuring non-negativity.]
  },
) <kmsz>

Steps 1--3 are over $QQ$ and produce rational $h_j$, $f$. The map $F_delta$ is non-negative by construction and non-SOS because $f$ contributes a non-SOS term, so the corresponding map is PNCP @Klep_2017.

=== Rationalization <rationalization>

Since step 4 of @kmsz involves a floating-point SDP, $delta$ is only known numerically. To obtain an exact rational certificate, we rationalize the Gram matrix post hoc.

Because $F_delta$ is non-SOS, we cannot directly find a Gram matrix for it. Instead, we relax: search for a Gram matrix $G$ for $F_delta dot S$, where $S = (sum_(i,j) x_i^2 y_j^2)^l$ is a fixed SOS multiplier. This bihomogeneous multiplier preserves the separate degrees in $bold(x)$ and $bold(y)$, matching the biquadratic structure of $F_delta$. Expanding $F_delta dot S$ in the monomial basis $bold(v)$ produces the feasibility problem
$ G succ.eq 0, quad chevron.l A_i, G chevron.r = b_i, quad i = 1, ..., m, $ <sos-sdp>
with rational data $A_i$, $b_i$ from the coefficient equations. If @sos-sdp is feasible, then $F_delta dot S$ is SOS, which certifies that $F_delta$ is non-negative @Klep_2017.

Let $G$ be a numerical solution with $mu = min "eig"(G) > 0$ and residual $epsilon = max_i |chevron.l A_i, G chevron.r - b_i|$. If $mu > epsilon$, a rational feasible $hat(G)$ is obtained by:
1. Exploiting the known $e$-dimensional nullspace @Klep_2017: set the corresponding blocks $tilde(G)_11 = 0$, $tilde(G)_12 = 0$.
2. Rationalize the remaining block $tilde(G)_22$ using continued fractions.
3. Project back to the affine constraint space $chevron.l A_i, G chevron.r = b_i$.
4. Verify $hat(G) succ.eq 0$ by a final rational feasibility SDP.

Once $hat(G)$ is rationalized, the coefficients of $F_delta$ are extracted by subtracting the contribution of the SOS multiplier $S$ from $bold(v)^T hat(G) bold(v)$. A final SDP check verifies that the extracted polynomial is not SOS, confirming the PNCP certificate is exact. Because IPMs produce solutions in the strict interior, the gap $mu > epsilon$ holds with room to spare, making the rationalization robust in practice.

= Methods <methods>
This chapter defines the computational workflow used in the current implementation. The goal is reproducibility with mathematical traceability: each computational step corresponds to a defined object or operation from the preceding chapter.

A note on normalization is in order. A quantum state is conventionally a PSD operator with unit trace, and a quantum channel is trace preserving. The tests in this thesis are, however, positivity- and sign-based: detecting entanglement amounts to checking the sign of a partial transpose, a witness expectation, or an SOS margin, none of which depend on the trace. We therefore work throughout with unnormalized operators, equivalently: with the convex cones of PSD and PPT matrices, and normalize only when a quantity must be interpreted as a physical state or channel. This avoids redundant rescaling in the inner loops without affecting any conclusion.

Two design constraints shape the workflow:
1. _exact_ separability testing is intractable in the dimensions where the conjecture is open (Section @dps-section), partial testing is very much doable, but every complete test is a one-sided relaxation;
2. SDP solutions are floating-point and fragile near feasibility boundaries (Section @rationalization), so positive detections must be re-validated exactly.

== Generating PPT candidates

We generate PPT states directly as the Choi matrices of PPT maps. A short sampling routine produces a random PPT state; an entanglement filter then keeps only the candidates that can plausibly yield a counterexample.

#algorithm-figure(
  "Generating a random PPT state",
  {
    Line[Sample a real matrix $R in M_(n m)$ with i.i.d. standard-normal entries and form the PSD matrix $rho = R R^T$.]
    Line[Compute the smallest eigenvalue of the partial transpose, $lambda = lambda_min (rho^Gamma)$. If $lambda < 0$, set $rho <- rho - lambda I$.]
  }
) <ppt-gen>

The construction is correct without any further work. By definition $rho = R R^T succ.eq 0$, and when $lambda < 0$ the shift adds the non-negative multiple $|lambda| I$ to $rho$, so positivity is preserved. The same shift acts on the partial transpose as $rho^Gamma - lambda I$, because $I^Gamma = I$, which raises every eigenvalue of $rho^Gamma$ by $|lambda|$ and hence makes $rho^Gamma succ.eq 0$. The resulting $rho$ is therefore PPT. Normal sampling is deliberate: an earlier prototype drew integer matrices for cosmetically nicer certificates, but those were not random enough and often had to be regenerated, whereas standard-normal entries almost always yield a valid candidate on the first try.

The off-diagonal blocks may optionally be symmetrized before the shift, for each block pair $(i, j)$ with $i < j$, both $rho_(i j)$ and $rho_(j i)$ are replaced by their average $(rho_(i j) + rho_(j i)) slash 2$, which forces $rho = rho^Gamma$, a fixed point of the partial transpose. This symmetrization is not essential; it only makes the computation slightly simpler and faster, since a single positivity shift then certifies both $rho succ.eq 0$ and $rho^Gamma succ.eq 0$ at once. Since the Choi matrix of a map is only defined up to the basis convention, restricting to such partial-transpose-invariant representatives is a permissible convenience rather than a requirement. The motivation is that our precomputed witnesses (Section @kmsz) are themselves generated in this partial-transpose-invariant shape, so matching the candidates to the same form may make them easier to detect; we therefore run the search in two versions, with the symmetry trick (the `ppt_invariant` option) and without it, in case the symmetric family is too narrow or escapes detection by the witnesses we have.

A random PPT state may be separable or bound entangled, and only the entangled ones are useful here. If either composed map is entanglement breaking, equivalently, has a separable Choi matrix, then the composition is automatically entanglement breaking and cannot violate the conjecture. We therefore discard separable candidates and keep only states for which we can certify entanglement, assembling a pool of genuine bound entangled states whose compositions are worth testing.

=== Alternative bound-entangled constructions

The random sampling above is fast, but it gives no control over entanglement: a sizable fraction of the samples are separable and thrown away by our filter, and the bound entangled survivors are unstructured. Several constructions in the literature instead produce bound entangled (PPT and entangled) states by design; we did not use them in the search, for the reasons noted below, but they are the natural starting points for a more structured candidate generation (cf. limitation 5 in @complexity-limits).

- *Unextendible product bases (UPB).* A UPB is a set of mutually orthogonal product vectors spanning a proper subspace whose orthogonal complement contains no product vector. The normalized projector onto that complement is PPT yet entangled @Bennett_1999. UPBs give explicit, low-rank bound entangled states, but the construction is the opposite of generic: it relies on specific, hand-picked bases that cannot be sampled at random, and in some dimensions may not exist at all, so it cannot drive a randomized search.
- *Antisymmetric-subspace states.* Sindici and Piani construct a simple class of PPT entangled states from the projectors onto the symmetric and antisymmetric subspaces of two identical systems, generalizing the Werner states @Sindici_2018. The construction is explicit, but certifying entanglement of the resulting states still reduces to an SDP, so it does not scale better than our pipeline.
- *Symmetric random induced states.* More recently, Louvet, Damanet, and Bastin study bound entanglement in symmetric random induced states @Louvet_2025. They sample a random pure state from the symmetric subspace $cal(H)_S^(N + N_A + 1)$ and trace out $N_A$ subsystems, leaving a mixed state on $cal(H)_S^(N + 1)$ that is, with high probability, bound entangled. The catch is dimensionality: for $N = 4$ bound entanglement is most likely at $N_A = 12$, which demands an enormous ancillary space and correspondingly heavy computation. Even then the probability of entanglement stays below $0.5$, comparable to the hit rate of our own construction. It is nonetheless a promising state-of-the-art source of candidates and a useful point of comparison for future work.

== Composing maps and testing for entanglement

Given two PPT candidates $Phi_1$, $Phi_2$ from the pool, we form the Choi matrix of their composition using the ampliation operation of Section @map-comp, $C_(Phi_1 compose Phi_2) = (I_n times.o Phi_1)(C_(Phi_2))$. We then test the composite for entanglement: if any composition is ever found to be entangled, the PPT2 conjecture is violated. Because composition is not commutative, the search ranges over _every ordered pair_ of pool states, self-pairs included.

Each composite is screened by three complementary detection tests:
- *Scalar witness test.* For each witness $W$ in the precomputed library, $"Tr"[W C] < 0$ certifies that $C$ is entangled.
- *Map witness test.* The _stronger_ condition $(I times.o Phi_W)(C) succ.eq.not 0$, evaluated as a smallest-eigenvalue check on the ampliation of the witness map; this can fire where the scalar test does not (Section @qse-section).
- *DPS relaxation.* The level-$2$ DPS relaxation (Section @dps-section); a feasible dual certificate flags entanglement.

#algorithm-figure(
  "Composition and entanglement test for one ordered pair",
  {
    Line[Take pool states $Phi_1, Phi_2$ and form $C = (I_n times.o Phi_1)(C_(Phi_2))$ by ampliation.]
    Line[Over the witness library, compute $min_W "Tr"[W C]$ and $min_W lambda_min ((I times.o Phi_W)(C))$.]
    Line[Compute the level-$2$ DPS robustness of $C$.]
    Line[Flag $C$ as entangled if any of the three scores indicates entanglement (beyond a tolerance), and report it for exact re-validation (Section @rationalization).]
  },
)

The witness library and the DPS test are complementary: precomputed witnesses are cheap to evaluate but are individually much less likely to detect entanglement, whereas the DPS relaxation is more robust but far more expensive. Running all three and recording each score also lets us compare the methods directly, how often the witnesses alone suffice versus how much the DPS relaxation adds (Section @results-section).

== Expanding the witness library <asym-witnesses>

The witness library of Stage 1 stores one matrix per polynomial, the symmetric Gram representative returned by the construction. As Section @gram-freedom shows, this is an arbitrary choice within the affine family $M_0 + L$, and every other representative is an equally valid witness over the reals. The earlier MATLAB prototype only ever used the symmetric representative; we instead exploit the freedom and expand each polynomial into several representations, enlarging the witness library at negligible cost. For each generated map $M_0$ we sample coefficients $lambda$ and form $M_0 + sum_alpha lambda_alpha N_alpha$ over a basis ${N_alpha}$ of $L$ (the `gram_freedom` primitive); the driver `gen_asym.jl` does this in bulk, emitting a configurable number of representations per source map.

The symmetric representative is also the partial-transpose-invariant one, $M_0 = M_0^Gamma$, the same shape as the symmetrized state construction of @ppt-gen. The added representatives break this symmetry. By construction they are unchanged as real polynomials, so they retain every property of the original witness on the real domain: they stay positive on all real product vectors and never assign a negative score to a separable state with real entries. They differ only in how they act on entangled inputs, which is precisely what makes them useful as additional, independent screens in the all-pairs search. A sampling-based block-positivity check (the `is_block_positive` primitive) guards the construction and is also the tool for probing the complex-domain question raised in Section @gram-freedom.

== Performance
A direct head-to-head comparison with the earlier MATLAB prototype is difficult, because the two implementations were built with different goals and make different design choices. The prototype drew random _integer_ matrices to produce cosmetically nicer certificates; this is not random enough and frequently failed, so a single candidate often needed many recomputation attempts. Our Julia implementation can still generate candidates that way, but for speed it samples from a normal distribution, which almost always succeeds on the first try (Section @ppt-gen). The two also rationalize differently: the prototype rationalizes the entire problem up front, paying for rational arithmetic throughout, whereas we rationalize only the final solution after the SDP is solved (Section @rationalization). The solver itself is largely language-independent, but the Julia code benefits from more efficient data handling and tighter solver integration.

With those caveats, we report rough wall-clock figures for the current implementation #margin-note[Add performance numbers]. For the same dimensions the prototype took orders of magnitude longer to produce candidates #margin-note[MATLAB perf computed/from article].

== Rationalization and Certificate Validation
Our constructions are intrinsically numerical, and the resulting maps are represented as floating-point matrices. This is sufficient for quick tests and exploration, but numerical noise can lead to false positives near the boundary of positivity. Therefore, we apply a post-solver rationalization and re-validation step to confirm any positive detections.

The rationalization, developed in Section @rationalization, proceeds as follows:
1. Extract the Gram matrix $G$ from the SDP constraints, $(delta f + h^2)(sum_(i,j) x_i^2 y_j^2)^l = bold(x)^T bold(y)^T G bold(y) bold(x)$.
2. Eigendecompose $G$, set the first $e = (n-1)(m-1)$ eigenvalues to zero, and rationalize the remaining coefficients to obtain $tilde(G)$.
3. Reformulating the constraints in terms of $tilde(G)$, isolate the coefficients of the polynomial $hat(p) = delta f + h^2$ from the relaxation terms $(sum_(i,j) x_i^2 y_j^2)^l$.
4. Solve a final feasibility SDP to verify that $hat(p)$ is _not_ a sum of squares, confirming the certificate is exact.

Upon completion we are confident in the validity of the certificate, so we can choose to represent it in floating point or as a rational number, as long as we are sure to use enough precision.

== The complete pipeline

The components above assemble into a three-stage pipeline. The first two stages build two pools independently, a library of PNCP witness maps and a pool of bound entangled PPT states, and the third stage tests all of their compositions. Separating generation from testing lets each pool be built once, checkpointed, and reused across runs.

#algorithm-figure(
  "PPT2 counterexample search",
  {
    Comment[Stage 1: witness library]
    For($i = 1, ..., N_W$, {
      Line[Construct a PNCP map via @kmsz, rationalize its certificate (Section @rationalization), and store the witness $W_i$.]
    })
    Comment[Stage 2: bound entangled candidate pool]
    For($j = 1, ..., N_C$, {
      Line[Sample a random PPT state via @ppt-gen; keep it only if we can verify it is entangled.]
    })
    Comment[Stage 3: test all compositions]
    For($"every ordered pair" (a, b)$, {
      Line[Form the composite channel $C = (I_n times.o Phi_a)(C_(Phi_b))$.]
      Line[Flag $C$ if any witness fires ($"Tr"[W_i C] < 0$ or $(I times.o Phi_(W_i))(C) succ.eq.not 0$) or the level-$2$ DPS relaxation reports entanglement.]
      Line[If $C$ is entangled, re-validate exactly and export the certificate: a counterexample to the PPT2 conjecture.]
    })
  },
) <pipeline>

In our experiments no composition ever produced an entangled composite; the search terminates after exhausting every ordered pair without a counterexample (@results-section).

== Implementation Architecture

The implementation is a small Julia package in the `code/` directory. The core logic lives in the `ppt2` module (`code/src/ppt2.jl`), with the PNCP construction split into `code/src/pncp.jl` and included into the same module. Together they map one-to-one onto the operations defined above:

#table(
  columns: 2,
  align: (left, left),
  stroke: none,
  table.header([*Function*], [*Role in the workflow*]),
  table.hline(),
  [`rand_ppt`], [Sample a random PPT state (@ppt-gen): the positivity shift, with optional block symmetrization (`ppt_invariant`).],
  [`ampliation`], [Compute $(I times.o Phi)(C)$, used both to compose maps and to apply a witness test.],
  [`sample_pncp_form`, `segre_kernel_basis`, `non_sos_form`], [The KMSZ construction (@kmsz): sample Segre-variety points, build the linear forms $h_j$, and produce the non-SOS quadratic form $f$.],
  [`solve_sos`], [Set up and solve the SOS feasibility/optimization SDP for a given relaxation degree $l$; optionally trigger rationalization.],
  [`rationalize_certificate`], [Post-solver rationalization (@rationalization): zero the first $e$ Gram eigenvalues, recover rational coefficients, and re-check non-SOS.],
  [`find_pncp_poly`, `pncp_mat`, `poly2mat`], [Orchestrate witness generation with retries and export the certificate as a Choi matrix.],
)

The module also exports a few supporting primitives: `rand_sep` and `rand_psd` for reference states, `is_ppt` for the PPT check, `gram_freedom` and `is_block_positive` for the witness-representation freedom and the block-positivity check of Section @gram-freedom, and `swap`/`antisymmetric_projector` for the antisymmetric-subspace construction of @Sindici_2018.

The package leans on the established Julia optimization and quantum-information stack rather than reimplementing it. Polynomials and the SOS cone are handled by `DynamicPolynomials` and `SumOfSquares`; the resulting semidefinite programs are modelled with `JuMP` and solved by `MOSEK` through `MosekTools` (any JuMP-compatible SDP solver could be substituted). The DPS hierarchy is not reimplemented: the search driver calls `entanglement_robustness` from `Ket`, an existing quantum-information toolbox that also supplies utilities such as the partial transpose. Matrices and witness libraries are serialized with `JLD2` so that generation and search can be separated and resumed.

Five command-line drivers in `code/scripts/` orchestrate the long-running jobs, all sharing a `common.jl` harness that provides resumable, reproducible, multithreaded batch generation: completed batches are detected and skipped on a rerun, and every candidate is seeded deterministically so a configuration yields the same dataset regardless of thread count. `gen_pncp.jl` builds the PNCP witness library (Stage 1), `gen_ppt.jl` samples and DPS-filters the bound entangled candidate pool (Stage 2), and `test_ppt2.jl` runs the threaded all-pairs search of @pipeline (Stage 3), logging any detection together with the offending state and witness. `gen_asym.jl` expands a witness library into the alternative Gram representations of Section @asym-witnesses. A fifth driver, `compare_detection.jl`, records every criterion's score on each detected state, the two witness tests and the DPS robustness, so the detection methods can be compared directly (Section @results-section). The `code/test/` suite checks the construction against reference values and verifies that generated maps are positive on large random samples, and the `code/notebooks/` directory documents the rationalization, PPT-state, and UPB workflows interactively.

== Complexity and Practical Limits <complexity-limits>
All the problems we are looking at are SDPs with exponential @Gharibian_2009 complexity in dimension and relaxation depth @Doherty_2004. Practical bottlenecks are:
1. processing power: these problems are by nature not easily parallelizable, simplex methods are intrinsically sequential and internal point methods rely on factoring large _sparse_ matrices, so the potential of GPU speedups is fairly low.
2. memory growth: the size of the SDP grows exponentially with dimension and relaxation depth, leading to memory bottlenecks even for moderate dimensions, i.e. for $4 times 4$ states DPS level 3 is already infeasible on standard hardware, requiring hundreds of gigabytes of RAM. Even improvements such as adding KKT constraints to achieve finite convergence @Harrow_2017 lead to significant increases in problem size, so they are not without cost.
3. solver instability near degeneracy: SDP solvers produce floating-point solutions, and when the problem is near the boundary of feasibility, small numerical errors can lead to incorrect conclusions about positivity. This is particularly problematic for our purposes, since we are interested in the very existence of such bound states.
4. Generated entanglement witnesses are not guaranteed to be optimal, there has been research on finding optimal witnesses @Lewenstein_2000, but this only applies to specific cases, so it remains open, how our witnesses can be improved. In fact many of our witnesses may even be decomposable, which would make them entirely unable to detect PPT-entangled states @Horodecki_2009.
5. Searching for PPT candidates in the first place is non-trivial, and random generation may not be sufficient to find counterexamples if they exist. We may require a more structured approach if the volume of the search space is in fact 0 (this would not mean the conjecture holds).

= Experimental Evaluation <results-section>
This section summarizes computational results from the pipeline implementations.

== Computational Results
All our tests were performed for the $4 times 4$ case, since it is the smallest dimension where the conjecture remains open. While our methods may handle even slightly higher dimensions, the computational cost grows rapidly, and we have not yet fully explored those regimes.

We have computed a library of 10,000 PNCP maps using the KMSZ construction. Each map was verified to be positive but not completely positive through rationalization and SOS checks.

Using Pipeline A (Randomized PPT Composition + DPS), we performed broad-coverage searches on random PPT-map compositions:
- 10,000 random candidates generated with PPT constraints
- DPS level 2 in addition to the precomputed entanglement witnesses
- Computational time: approximately 5 hours on 80 CPU cores for complete verification

No counterexamples to the PPT2 conjecture were found in these searches.

== Representation and symmetry effects <witness-findings>

Alongside the main search we recorded, for every detected state, the score of each detection test, which lets us compare the methods directly. Two robust patterns emerge.

First, the map witness test is empirically the stronger one, exactly as argued in Section @qse-section. Across all detected composites, every state flagged by the scalar test $"Tr"[W C] < 0$ is also flagged by the map test $lambda_min ((I times.o Phi_W)(C)) < 0$, while the converse fails: the map test detects a strict superset. The two run together and their scores are logged by `compare_detection.jl`.

Second, the partial-transpose-invariant (symmetric) representatives are inert for our purpose. Symmetrized witnesses, the $M_0 = M_0^Gamma$ representatives used by the original construction, never detected a single PPT entangled state in our experiments, and PPT states built in the symmetric `ppt_invariant` shape never appeared as counterexamples. This is consistent with the fact that decomposable witnesses cannot detect PPT entanglement (Section @qse-section): a witness aligned with the partial-transpose fixed-point set carries little information about the PPT cone. We do not have a full characterization, but the observation is uniform enough that we run the search on the asymmetric representatives (Section @asym-witnesses) and on generic, non-symmetrized candidates instead.

The asymmetric representatives behave exactly as the theory predicts on the real domain. They remain positive on all real product vectors, and over many random trials they never assigned a negative score to a real separable state, so they raise no false positives in the regime we actually search, where every generated state and witness has real entries, and they detect real-domain entanglement readily. Whether this extends to the complex domain is open: the relations underlying the representation freedom (Section @gram-freedom) vanish only on real product vectors, so an asymmetric representative that is a valid witness over $RR$ need not be block-positive over $CC$. Establishing complex block-positivity, or restoring it by a suitable projection, is left for future work.

= Conclusion
This thesis presents a computational approach to the PPT2 conjecture, integrating methods for generating PPT maps, composing them, and testing for entanglement using both DPS-based relaxations and PNCP-based witness construction. While no counterexamples were found in our extensive searches, the methodology provides a framework for further exploration and potential discovery in this open problem.

#bibliography(
  title: "References",
  style: "ieee",
  "bibliography.bib"
)
