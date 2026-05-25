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
  abstract_en: [This thesis explores...],
  abstract_sl: [V tem delu raziskujemo...],
  extended_abstract_sl: [Daljši slovenski povzetek vsebine...],
)

#set heading(numbering: "1.1")
#show heading: it =>[
    #block(it.body)
]
//#show heading.where(level: 1): it =>[
//    #block(it.body)
//]

#set math.equation(numbering: "(1)")

#let (
  theorem, remark, example, definition, proof, rules: thm-rules
) = default-theorems("thm-group", lang: "en", thm-numbering: thm-numbering-linear)
#show: thm-rules

#let reference = margin-note("reference")
#let PPT2 = $"PPT"^2$

= TODO <todo>
- [ ] write introduction - what is entanglement, how do we use it, what is the PPT2 conjecture and why it matters, what are its implications if true or false. Current state of knowledge i.e. known results for $n<=3$ expectations for $n>=4$, testing entanglement is hard...
- [ ] generating PPT candidates, 1. start with psd matrices $A A^T$ (optionally symmetrize) and add smallest eigenvalue of partial transpose on diagonal to ensure PPT - the downside here is that the state may not be entangled, 2. use unextendible product bases (UPB) - this should guarantee bound entanglement, i.e. PPT and entangled, requires references to UPB papers and some explanation
- [ ] A note on composition calculation; not as straightforward as multiplying matrices, needs the Choi-Jamiolkowski transformation so that we can work with matrices, define _ampliation_ $(I times.o C_Phi)(C_Psi)$ which we use to compute the Choi matrix of the composition, as well check entanglement 
- [ ] DPS-based entanglement detection, describe hierarchy, SDPs in general, witness extraction for comparison later, extensions of DPS, i.e. adding KKT constraints, our code uses existing implementations (cite Ket library)
- [ ] Abhisheks PnCP map generation method, describe in more detail what it does and how it works, Julia implementation does rationalization differently for performance reasons, i.e. only after computing SDP, explain why this still works, get matlab and julia performance numbers to compare, explain how the matrix representation is chosen for each map (not unique), chosen to be invariant to partial transpose, but doesn't need to be - discussion
- [ ] get references for all statements, known results, methods, etc.
- [ ] write out entire algorithm/pipeline - loop: 1. generate candidate(s), 2. compute composition, 3. test for entanglement (using DPS and/or precomputed witnesses), 4. if candidate is entangled, rationalize and retest to confirm, 5. if confirmed, export certificate and stop
- [ ] results: no counterexamples found on 10000 random candidates, in 4x4 this takes about five hours to verify using DPS level 2 and a collection of 10000 precomputed witnesses
- [ ] issues with generated entanglement witnesses: almost surely not optimal, possibly even decomposable, so probably fails to detect most entangled states
- [ ] issues with DPS: numerical instability, only lower levels of hierarchy are computationally feasible
- [ ] required sections: 1. introduction - explains entanglement, PPT2 conjecture 2. theoretical background - required formalism: positivity, complete positivity, entanglement, partial transpose, quantum states/channels, entanglement witnesses, DPS hierarchy, semidefinite programming, Choi representation, maps/matrices/operators 3. methods - generating candidates, why we can't use known constructions, compositions, testing entanglement via DPS and witnesses, generating witnesses via PNCP maps, post-solver numerical validation, performance notes 4. results/discussion/conclusion - TBD

#inline-note[careful with "improvements" on the DPS hierarchy; some of them work only for specific types of states i.e. diagonal unitary invariant states, for which the #PPT2 conjecture has already been proven]

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
This chapter establishes the mathematical framework used throughout the thesis. We fix notation, define the key objects, and develop the tools needed for the computational pipeline in Chapter 3.

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

Every CP map is positive. Maps that are positive but not completely positive (PNCP) are central to entanglement detection.

#definition(name: "Choi-Jamiolkowski isomorphism")[
  The _Choi matrix_ of a linear map $Phi: M_n -> M_m$ is
  $ C_Phi = (I_n times.o Phi)(omega_n) = sum_(i,j=1)^n E_(i j) times.o Phi(E_(i j)) in M_n times.o M_m, $
  where $omega_n = sum_(i,j) E_(i j) times.o E_(i j)$ is the (unnormalized) maximally entangled state. The assignment $Phi arrow.bar C_Phi$ is a linear isomorphism between maps $M_n -> M_m$ and matrices in $M_n times.o M_m$.
]
#inline-note[the definition of $omega$ here is unnecessary. Reformulate the the isomorphism for more clarity. It is important that both the map and the matrix forms are clearly stated though.]

Under this isomorphism, $Phi$ is CP if and only if $C_Phi succ.eq 0$ @Choi_1975, and $Phi$ is trace-preserving if and only if $"Tr"_B [C_Phi] = I_n$.

#definition(name: "Partial transpose")[
  For $rho = (rho_(i j))_(i,j=1)^m in M_m times.o M_n$ in block form with $rho_(i j) in M_n$, the _partial transpose_ with respect to subsystem $B$ is
  $ rho^Gamma := (I_m times.o T)(rho) = (rho_(i j)^T)_(i,j=1)^m, $
  where $T$ is the transposition map. For Hermitian $rho$, the choice of subsystem does not affect whether $rho^Gamma succ.eq 0$, so we drop the subscript.
]
#inline-note[This definitions is a bit clunky, Instead explicitly state $rho$ as a kronecker product of two matrices. That way it is clear what we mean by 'subsystem'. The final sentence is also a bit confusing, rephrase to clarify that for Hermitian matrices, the partial transpose with respect to either subsystem is the same, so we can just write $rho^Gamma$ without specifying the subsystem.]

#definition(name: "PPT state and PPT map")[
  A state $rho in M_m times.o M_n$ is _PPT_ if $rho succ.eq 0$ and $rho^Gamma succ.eq 0$. A map $Phi: M_n -> M_m$ is _PPT_ if its Choi matrix is PPT:
  $ C_Phi succ.eq 0 quad "and" quad C_Phi^Gamma succ.eq 0. $
]

Since $C_Phi succ.eq 0$ characterizes CP maps, a PPT map is a CP map whose Choi matrix has non-negative partial transpose. The composition of two PPT maps is trivially PPT.

=== Map composition via Choi matrices

The Choi matrix of a composition $Phi compose Psi$ is not the product $C_Phi C_Psi$. Expanding the definition directly gives:
$
  C_(Phi compose Psi)
  = sum_(i,j) E_(i j) times.o Phi(Psi(E_(i j)))
  = sum_(i,j) E_(i j) times.o Phi(sum_(k,l) (C_Psi)_(i k, j l) E_(k l)) \
  = sum_(i,j,k,l,p,q) (C_Psi)_(i k, j l) (C_Phi)_(k p, l q) E_(i j) times.o E_(p q),
$
yielding the index-wise formula
$ (C_(Phi compose Psi))_(i p, j q) = sum_(k,l) (C_Psi)_(i k, j l) (C_Phi)_(k p, l q). $
#inline-note[This is nothing new just a note on the pitfalls of maps and matrices. Simplify the notation a bit and add in the ampliation operation to relate it to other operations we use later on.]

== Quantum States and Entanglement

A _quantum state_ on $CC^n$ is a density matrix: a PSD matrix $rho in M_n^+$ with $"Tr"[rho] = 1$. A _quantum channel_ is a CPTP map $Phi: M_n -> M_m$.

#inline-note[state what CPTP means: completely positive trace preserving]

#definition(name: "Separability and entanglement")[
  A bipartite state $rho in M_m times.o M_n$ is _separable_ if
  $ rho = sum_i p_i rho_i^A times.o rho_i^B, quad p_i >= 0, quad sum_i p_i = 1, quad rho_i^A in M_m^+, quad rho_i^B in M_n^+. $
  Otherwise $rho$ is _entangled_.
]

Deciding separability is NP-hard in general @Gharibian_2009, which motivates necessary conditions and convex relaxations.

#definition(name: "Entanglement-breaking map")[
  A map $Phi: M_n -> M_m$ is _entanglement breaking_ (EB) if $(I_k times.o Phi)(rho)$ is separable for every $k in NN$ and every state $rho in M_k times.o M_n$.
]
#inline-note[Find a reference for this definition. Is the use of $k$ and ampliation necessary?]

=== Separability criteria

Several necessary conditions for separability allow entanglement detection when violated. #inline-note[what does this mean?]

#inline-note[the article "Entanglement detection" by O. Guhne and G. Toth https://doi.org/10.1016/j.physrep.2009.02.004 provides a good overiview of separability criteria. Use it as a reference and source for the following paragraphs, and add more details and references as needed.]

*PPT criterion.* If $rho$ is separable, then $rho^Gamma succ.eq 0$ @Peres_1996 @Horodecki_1996. The condition is necessary but not sufficient: bound entangled (PPT entangled) states exist whenever $m n > 6$ @Horodecki_2009.

The PPT criterion is an example of a criterion using a positive but not completely positive map. #inline-note[this is an important statement to make, but put it in context.]

*Range criterion.* If $rho in M_m times.o M_n$ is separable, there exist product vectors $lr({|psi_i chevron.r times.o |phi_i chevron.r})$ spanning the range of $rho$ such that $lr({|psi_i chevron.r times.o overline(|phi_i chevron.r)})$ spans the range of $rho^Gamma$ @Horodecki_1997. The range criterion can detect certain PPT entangled states, but fails when $rho$ is full rank (e.g., under noise), since then any set of vectors spans its range trivially.

*Reduction criterion.* If $rho in M_m times.o M_n$ is separable, then
$ rho_A times.o I_n - rho succ.eq 0 quad "and" quad I_m times.o rho_B - rho succ.eq 0, $
where $rho_A = "Tr"_B [rho]$ and $rho_B = "Tr"_A [rho]$ are the marginals @Horodecki_1999. #inline-note[If we use the word marginals, we must also define it.]

#inline-note[CCNR criterion]

*Majorization criterion.* If $rho$ is separable, then $lambda(rho) prec lambda(rho_A)$ and $lambda(rho) prec lambda(rho_B)$, where $lambda(dot)$ denotes the non-increasingly ordered eigenvalue vector and $a prec b$ denotes majorization @Nielsen_2001. This follows from the reduction criterion @Nielsen_2001 and shares its limitations. #inline-note[clarify what majorization is, and use a different symbol for it, since $prec$ and $succ$ are already used for positive definiteness.]

#inline-note[Additional criteria of interest include:]
Algorithmic approaches reformulate separability as a convex optimization problem or a semidefinite program. One of the most powerful separability criteria is the method of symmetric extensions.

Criteria based on covariance matrices can be very strong when they work, but it is not clear how to use them in the multipartite case.

Separability criteria have also been developed using linear contractions and permutations, leading to the PPT and CCNR criteria.

Other approaches include the use of special observables, either bell inequalities or entanglement witnesses.

=== Entanglement witnesses

#definition(name: "Entanglement witness")[
  A Hermitian operator $W$ is an _entanglement witness_ if $"Tr"[W sigma] >= 0$ for all separable $sigma$, and $"Tr"[W rho] < 0$ for some entangled $rho$.
]

By the Hahn-Banach theorem, for every entangled state $rho$ there exists a witness detecting it, so the family of all witnesses completely characterizes the set of separable states @Horodecki_2009. Witnesses are one-sided: a negative expectation value certifies entanglement, but a non-negative one does not certify separability.

#inline-note[what is the Hahn-Banach theorem?]

#definition(name: "Decomposable witness")[
  A witness $W$ is _decomposable_ if $W = P + Q^Gamma$ for some $P, Q succ.eq 0$; otherwise it is _non-decomposable_ @Lewenstein_2000.
]

Decomposable witnesses cannot detect PPT entangled states: for any PPT state $sigma$, $"Tr"[(P + Q^Gamma) sigma] = "Tr"[P sigma] + "Tr"[Q sigma^Gamma] >= 0$. Only non-decomposable witnesses are useful for our purposes.

#inline-note[add figure showing geometric representation of entanglement witnesses]

Under the Choi-Jamiolkowski isomorphism, every entanglement witness $W = C_Phi$ corresponds to a PNCP map $Phi$. Note that the map condition $(I_k times.o Phi)(rho) succ.eq.not 0$ is strictly stronger than the scalar condition $"Tr"[C_Phi rho] < 0$ @Horodecki_2009.

== Semidefinite Programming

A _semidefinite program_ (SDP) is a convex optimization problem in which a linear objective is minimized subject to a linear matrix inequality (LMI) constraint. The standard primal form is

$ "minimize" &quad c^T bold(x) \
  "subject to" &quad F(bold(x)) := F_0 + sum_(i=1)^n x_i F_i succ.eq 0, $ <sdp-primal>

where $bold(x) in RR^n$ is the optimization variable, $c in RR^n$ is the cost vector, and $F_0, F_1, ..., F_n in H_d$ are fixed Hermitian matrices. The associated _dual SDP_ is

$ "maximize" &quad -"Tr"[F_0 Z] \
  "subject to" &quad Z succ.eq 0, \
               &quad "Tr"[F_i Z] = c_i, quad i = 1, ..., n. $ <sdp-dual>

When the Slater #margin-note[what is the Slater condition? Either define or reference it.] condition holds (strictly feasible primal or dual), strong duality holds and primal and dual optima coincide. When $c = 0$, @sdp-primal is a _feasibility problem_: determine whether there exists $bold(x)$ with $F(bold(x)) succ.eq 0$. If infeasible, a dual feasible $Z succ.eq 0$ satisfying $"Tr"[F_i Z] = 0$ with $"Tr"[F_0 Z] > 0$ certifies infeasibility.

=== Interior point methods

#inline-note[we are only talking about interior point methods because the solvers we (and most state-of-the-art solvers) use are base on them. There are also other methods. Specify any quirks of this particular method that are relevant for our implementation. Reference MOSEK for details on the solver we use.]
SDPs are solved in practice by _interior point methods_ (IPMs), which follow a central path through the strict interior of the feasible region. An IPM applied to a $d times d$ SDP with $n$ variables reaches $epsilon$-accuracy in $O(sqrt(d) log(1 slash epsilon))$ iterations, each costing $O(n^2 d^2 + n d^3)$ arithmetic operations @Vandenberghe_1996. Because solvers produce floating-point solutions, feasibility conclusions near the boundary require post-solver validation (see Section @rationalization).

=== Entanglement testing as an SDP

Testing whether a state $rho$ is separable can be formulated as
$ h_("SEP")(M) = max {"Tr"(M rho) : rho in "SEP"}, $
which is NP-hard @Harrow_2017. #margin-note[wrong reference] The DPS hierarchy replaces SEP with a sequence of SDP-representable outer approximations that converge to it, giving a tractable relaxation at each level.

#inline-note[@Harrow_2017 actually provides a different relaxation $h_("ProdSym")$, which can be mapped to $h_("SEP")$.]

== DPS Hierarchy

The Doherty-Parrilo-Spedalieri (DPS) hierarchy @Doherty_2004 provides a complete sequence of SDP relaxations of separability. #inline-note[complete in the sense that for any entangled state, there exists a level of the hierarchy that detects it. In theory, with enough computational resources, we could test all levels and detect all entangled states. In practice, only the first few levels are computationally feasible, so we can only detect a subset of entangled states.]

=== Symmetric extensions

For a separable state $rho = sum_i lambda_i x_i x_i^* times.o y_i y_i^*$, define the $k$-fold _symmetric extension_
$ rho_k = sum_i lambda_i x_i x_i^* times.o (y_i y_i^*)^(times.o k) in H(CC^m times.o (CC^n)^(times.o k)). $
This extension satisfies
$ (I_m times.o Pi_k) rho_k (I_m times.o Pi_k) = rho_k, quad rho_k^(Gamma_s) succ.eq 0 quad forall s = 1,...,k, quad "Tr"_(B_2 ... B_k)[rho_k] = rho, $
where $Pi_k$ is the projector onto the bosonic (fully symmetric) subspace of $(CC^n)^(times.o k)$.
#inline-note[this is somewhat technical. Maybe state each condition separately and explain its meaning.]

#definition(name: "DPS set at level k")[
  $ "DPS"_n^k = lr({rho in H(CC^m times.o CC^n) : exists rho_k in H(CC^m times.o (CC^n)^(times.o k)) "satisfying the above three conditions"}). $
] #inline-note[instead of "above three conditions", refer to the specific equations]

Each $"DPS"_n^k$ is defined by SDP constraints, so membership is testable in polynomial time for fixed $k$. The hierarchy satisfies @Doherty_2004:
1. $S E P_n subset.eq "DPS"_n^k$ and $"DPS"_n^(k+1) subset.eq "DPS"_n^k$ for all $k >= 1$.
2. $"DPS"_n^1$ is equivalent to the PPT criterion.
3. Asymptotic completeness: $inter.big_(k >= 1) "DPS"_n^k = S E P_n$.

=== DPS as an explicit SDP

Testing $rho in "DPS"_n^k$ is the following feasibility SDP in the variable $rho_k$:
$ "find" &quad rho_k succ.eq 0 \
  "s.t." &quad (I_m times.o Pi_k) rho_k (I_m times.o Pi_k) = rho_k, \
         &quad rho_k^(Gamma_s) succ.eq 0, quad s = 1, ..., k, \
         &quad "Tr"_(B_2 ... B_k) [rho_k] = rho. $

#inline-note[Avoid re-stating the same formulas over and over again.]

If feasible, $rho$ passes the level-$k$ test. If infeasible, the dual certificate $Z succ.eq 0$ yields an entanglement witness for $rho$ via $W = -"Tr"_(B_2 ... B_k) [Z]$.

#inline-note[This part is important, we would like to obtain a witness, as that helps verify our result is not just a numerical fluke. Check the article and add more details on how the witness is obtained from the dual solution.]

=== Improvements and limitations

Adding KKT conditions to the SDP achieves _finite convergence_ at some level $k$ @Harrow_2017, but increases problem size and complicates witness extraction. Specialized hierarchies for diagonal unitary invariant states @Britz_2025 cannot be used here since the PPT2 conjecture is already proven for that family @Singh_2022.

#inline-note[Expand on this paragraph, be more specific about what each improvement is and how it works. Then it will be more clear why we didn't use them. Also search for any other improvements that might be relevant.]

== Sum-of-Squares Polynomials and PNCP Maps

A complementary approach to DPS exploits the polynomial representation of linear maps and the connection between sum-of-squares relaxations and SDPs.

#inline-note[what is the connection between SOS and SDPs?]

=== Non-negative polynomials and SOS

Let $RR[bold(x), bold(y)]$ be the ring of real polynomials in $bold(x) in RR^n$, $bold(y) in RR^m$. A polynomial $p$ is _non-negative_ if $p(bold(x), bold(y)) >= 0$ for all inputs; it is a _sum of squares_ (SOS) if $p = sum_i q_i^2$ for polynomials $q_i$. Every SOS polynomial is non-negative, but not every non-negative polynomial is SOS — the gap between the two gives rise to PNCP maps.

=== Gram matrix and SDP connection

#definition(name: "Gram matrix representation")[
  A homogeneous polynomial $p$ of degree $2d$ is SOS if and only if there exists $G succ.eq 0$ such that
  $ p(bold(z)) = bold(v)(bold(z))^T G bold(v)(bold(z)), $
  where $bold(v)(bold(z))$ is the vector of monomials of degree $d$ in $bold(z)$.
]

Testing whether $p$ is SOS therefore reduces to a semidefinite feasibility problem: find $G succ.eq 0$ satisfying the linear constraints that equate the coefficients of $bold(v)^T G bold(v)$ with those of $p$. If no such $G$ exists, $p$ is not SOS. #margin-note[do not reiterate the formula, instead refer to the previous definition.]

=== Polynomial representation of maps

Each linear map $Phi: M_n -> M_m$ corresponds to a biquadratic polynomial:
$ p_Phi (bold(x), bold(y)) := bold(y)^T Phi(bold(x) bold(x)^T) bold(y). $

The fundamental correspondence @Klep_2017 is:
- $Phi$ is _positive_ if and only if $p_Phi >= 0$ on $RR^n times RR^m$.
- $Phi$ is _completely positive_ if and only if $p_Phi$ is SOS.

#inline-note[is non-negative the same as positive? Be consistent with our terminology.]

A PNCP map therefore corresponds exactly to a non-negative non-SOS polynomial. This gives a third equivalent representation of maps (alongside matrices and linear operators) and allows entanglement witnesses to be constructed via polynomial optimization.

#algorithm-figure(
  "KMSZ construction for PNCP maps",
  {
    Line[Generate random points $x^((1)), ..., x^((t)) in RR^n$ and $y^((1)), ..., y^((t)) in RR^m$.]
    Line[Form bilinear forms $h_j (bold(x), bold(y)) = chevron.l x^((j)), bold(x) chevron.r dot chevron.l y^((j)), bold(y) chevron.r$ for $j = 1, ..., t$.]
    Line[Find $f in.not "span"{h_1, ..., h_t}$, and not SOS.]
    Line[Solve the SDP: find the smallest $delta > 0$ such that $F_delta = delta f + sum_j h_j^2 >= 0$.]
  },
) <kmsz>

#inline-note[this is a very high level description of the algorithm, add some more details on what each step means, but reference the original article for the technical details.]

Steps 1--3 involve only linear algebra over $QQ$ and produce rational $h_j$, $f$. The map $F_delta$ is non-negative by construction and non-SOS because $f$ contributes a non-SOS term, so the corresponding map is PNCP @Klep_2017.

=== Rationalization <rationalization>

Since step 4 of @kmsz solves a floating-point SDP, $delta$ is only known numerically. To obtain an exact rational certificate we rationalize the Gram matrix post hoc. The SDP arising from the SOS relaxation is the feasibility problem
$ G succ.eq 0, quad chevron.l A_i, G chevron.r = b_i, quad i = 1, ..., m, $
with rational data $A_i$, $b_i$.

#inline-note[clarify how we obtain this SDP from the SOS relaxation. Either write out the SOS relaxation before, or reference the relevant article where this is explained in detail.]

Let $G$ be a numerical solution with $mu = min "eig"(G) > 0$ and residual $epsilon = max_i |chevron.l A_i, G chevron.r - b_i|$. If $mu > epsilon$, a rational feasible $hat(G)$ is obtained by:
1. Exploiting the known $e$-dimensional nullspace @Klep_2017: set the corresponding blocks $tilde(G)_11 = 0$, $tilde(G)_12 = 0$.
2. Rationalize the remaining block $tilde(G)_22$ using continued fractions.
3. Project back to the affine constraint space $chevron.l A_i, G chevron.r = b_i$.
4. Verify $hat(G) succ.eq 0$ by a final rational feasibility SDP.

#inline-note[This rationalization procedure is taken from the article "there are many more positive maps than completely positive maps" @Klep_2017, but what we do is a little different. We Rationalize $G$ in the same way. However, we are generating a polynomial $F_delta$ that is not a sum of squares. Verifying it is positive is therefore hard (reference). Instead we relax the problem and verify that $F_delta$ times some SOS relaxation term is a sum of squares. (ref why this works) Once we obtain the rationalized $hat(G)$, we extract the terms of $F_delta$ by removing the SOS relaxation terms. So, the final feasibility SDP is to verify the resulting polynomial is not a sum of squares.]

#inline-note[The final feasibility SDP is only to verify that the rationalization process did not introduce any errors. The SDP is not "rational" in the sense that it does not have rational data. However, since IPMs yield results strictly in the interior of the feasible region, we can be confident in the results.]

== PPT2 Conjecture

#theorem(name: "PPT2 Conjecture")[
  If $Phi_1$ and $Phi_2$ are PPT maps, then $Phi_1 compose Phi_2$ is entanglement breaking @Christandl_2019.
]<ppt2>

The conjecture is proven for $n = 2$ (trivially, since all PPT states are separable in that dimension #margin-note[reference or explanation]), for $n = 3$ @Chen_2019, and for Choi-type maps in all dimensions @Singh_2022. It remains open for $n >= 4$ in the general case.

#inline-note[move this section up as far as possible, so that it is clear from the start what we are trying to do and why. We need to define linear maps, PPT and such beforehand.]

= Methods <methods>
This chapter defines the computational workflow used in the current implementation.

#inline-note[Although channels are trace preserving by definition, the computational tests in this thesis are primarily positivity/sign based. Therefore, intermediate computations may use non-normalized maps, then normalize where needed for interpretation.]

#inline-note[A quantum state is usually represented by a positive semidefinite operator normalized with unit trace. In this work, we mainly work with unnormalized quantum states and consider its convex cone.]

The goal is reproducibility with mathematical traceability: each code step corresponds to a defined object or operation from Chapter 2.

Design constraints:
1. separability testing difficulty in relevant dimensions,
2. numerical fragility near SDP feasibility boundaries,

== PPT Maps

#algorithm-figure(
  "Generating PPT maps",
  {
    Line[Generate a matrix $R$ with entries sampled from a normal distribution.]
    Line[$A = R^T R$ is a positive semidefinite matrix.]
    Line[Optionally symmetrize each block of $A$ to ensure the resulting map is invariant under partial transpose.]
    Line[Adjust the spectrum of $A$ to ensure the resulting map is PPT: $A_(P P T) = A - lambda I$, where $lambda < 0$ is the smallest eigenvalue of $A^Gamma$.]
  }
)

== Composite maps and entanglement testing
Compute Composite PPT maps from generated candidates and test with DPS and precomputed witnesses

#algorithm-figure(
  "Randomized PPT Composition + DPS",
  {
    Line[Generate random PPT maps $Phi_1$, $Phi_2$ using the method described above.]
    Line[Compute the Choi matrix of the composition $C_(Phi_1 compose Phi_2)$ using the ampliation operation.]
    Line[Test the resulting state for entanglement using DPS hierarchy and/or precomputed witnesses.]
  },
)

== Performance
Our performance gains stem from many factors. The MATLAB implementation originally generated random integer matrices to produce nicer results, but the downside is that this aproach is not random enough, and often failed. Resulting in a single candidate requiring many recomputation attempts. The Julia implementation allws for such candidate generation, but with speed in mind we opted for the more reliable approach of sampling from a normal distribution, which is far more likely to produce a valid candidate on the first try. While the solver's performance is generally independent of the language, the Julia implementation benefits from more efficient data handling and better integration with the solver, leading to faster overall runtimes. Additionally the rationalization steps in MATLAB are noticeably more expensive, both due to the overhead of rational arithmetic and the choice to rationalize the entire problem beforehand, while in Julia we only rationalize the solution after solving the SDP. All this to say, direct comparison between implementations is somewhat difficult and arguably pointless. It would be like comapring apples to oranges. With that in mind we present some rough performance numbers for the current implementation: #margin-note[Add performance numbers], The previous implementation was made with other goals in mind, but producing candidates for the same dimensions took orders of magnitude more time: #margin-note[MATLAB perf computed/from article].

== Rationalization and Certificate Validation
Our constructions are intrinsically numerical, and the resulting maps are represented as floating-point matrices. This is sufficient for quick tests and exploration, but numerical noise can lead to false positives near the boundary of positivity. Therefore, we apply a post-solver rationalization and re-validation step to confirm any positive detections.

The rationalization as explained previously is implemented in the `_fix_gram` function and works as follows:
1. Extract Gram matrix $G$ from constraints of the SDP $(delta f + h^2)(sum x_i^2 y_j^2)^l = bold(x)^T bold(y)^T G bold(y) bold(x)$.
2. Compute eigendecomposition of $G$, set the first $e$ eigenvalues equal to zero, and rationalize the remaining coefficients to obtain $tilde(G)$.
3. Reformulating the constraints in terms of $tilde(G)$, isolate the coefficients of polynomial $hat(p) = delta f + h^2$ from the relaxation terms $(sum x_i^2 y_j^2)^l$.
4. Solve another feasibility SDP to verify $hat(p)$ is not a sum of squares.

== Implementation Architecture
The active implementation is in the `code` directory, centered on `code/src/ppt2.jl`. #inline-note[details on code architecture - main stuff in src folder, additional scripts for reproducibility, experiments notebooks, dependencie - general julia stuff, any JuMP compatible solvers, Ket as existing quantum toolbox, and usage instructions]

== Complexity and Practical Limits
All the problems we are looking at are SDPs with exponential @Gharibian_2009 complexity in dimension and relaxation depth @Doherty_2004. Practical bottlenecks are:
1. processing power: these problems are by nature not easily parallelizable, simplex methods are intrinsically sequential and internal point methods rely on factoring large _sparse_ matrices, so the potential of GPU speedups is fairly low. 
2. memory growth: the size of the SDP grows exponentially with dimension and relaxation depth, leading to memory bottlenecks even for moderate dimensions, i.e. for $4 times 4$ states DPS level 3 is already infeasible on standard hardware, requiring hundreds of gigabytes of RAM. Even improvements such as adding KKT constraints to achieve finite convergence @Harrow_2017 lead to significant increases in problem size, so they are not without cost.
3. solver instability near degeneracy: SDP solvers produce floating-point solutions, and when the problem is near the boundary of feasibility, small numerical errors can lead to incorrect conclusions about positivity. This is particularly problematic for our purposes, since we are interested in the very existence of such bound states.
4. Generated entanglement witnesses are not guaranteed to be optimal, there has been research on finding optimal witnesses @Lewenstein_2000, but this only applies to specific cases, so it remains open, how our witnesses can be improved. In fact many of our witnesses may even be decomposable, which would make them unable to detect PPT-entangled states @Horodecki_2009.
5. Searching for PPT candidates in the first place is non-trivial, and random generation may not be sufficient to find counterexamples if they exist. We may require a more structured approach if the volume of the search space is in fact 0 (this would not mean the conjecture holds).

= Experimental Evaluation
This section summarizes computational results from the pipeline implementations.

== Computational Results
All our tests were performed for the $4 times 4$ case, since it is the smallest dimension where the conjecture remains open. While our methods may handle even slightly higher dimensions, the computational cost grows rapidly, and we have not yet fully explored those regimes.

We have computed a library of 10,000 PNCP maps using the KMSZ construction. Each map was verified to be positive but not completely positive through rationalization and SOS checks.

Using Pipeline A (Randomized PPT Composition + DPS), we performed broad-coverage searches on random PPT-map compositions:
- 10,000 random candidates generated with PPT constraints
- DPS level 2 in addition to the precomputed entanglement witnesses
- Computational time: approximately 5 hours on 80 CPU cores for complete verification

No counterexamples to the PPT2 conjecture were found in these searches.

= Conclusion
This thesis presents a computational approach to the PPT2 conjecture, integrating methods for generating PPT maps, composing them, and testing for entanglement using both DPS-based relaxations and PNCP-based witness construction. While no counterexamples were found in our extensive searches, the methodology provides a framework for further exploration and potential discovery in this open problem.

= Links
- https://arxiv.org/pdf/1506.08834
- https://arxiv.org/pdf/quant-ph/0603199
- https://arxiv.org/pdf/2402.12944
- https://arxiv.org/pdf/1807.01266
- https://arxiv.org/pdf/2001.01181
- https://arxiv.org/pdf/2010.07898
- https://arxiv.org/pdf/1807.03636
- https://arxiv.org/pdf/2011.03809
- https://arxiv.org/pdf/2512.06551
- https://arxiv.org/html/2501.03959v1
- https://arxiv.org/pdf/quant-ph/0308032
- https://felixleditzky.info/teaching/ST23/Felix%20Leditzky%20-%20Math%20595%20Quantum%20channels.pdf
- https://arxiv.org/pdf/quant-ph/0602223
- https://arxiv.org/pdf/0810.4507
- https://arxiv.org/pdf/1309.7992
- https://journals.aps.org/pra/abstract/10.1103/PhysRevA.71.032333
- https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.023101
- https://arxiv.org/pdf/quant-ph/9801069
- https://arxiv.org/pdf/quant-ph/0702225
- https://arxiv.org/pdf/0907.4979
- https://link.springer.com/article/10.1007/s00220-017-2859-0
- http://congres.cran.univ-lorraine.fr/2002/CDC_2002/pdffiles/papers/832_FrP06-6.pdf
- https://arxiv.org/html/2506.11346v3
- https://arxiv.org/html/2308.07019v5#bib.bib22
- https://helda.helsinki.fi/server/api/core/bitstreams/486c1a74-82f8-4c6b-b665-6a9142606c78/content
- https://arxiv.org/pdf/0712.1114
- https://www.nature.com/articles/s41598-022-14920-5
- https://arxiv.org/pdf/0907.2369
- https://iopscience.iop.org/article/10.1088/0305-4470/39/45/020/pdf
- https://arxiv.org/pdf/0805.1318
- https://iopscience.iop.org/article/10.1088/1751-8121/acaa16
- https://link.springer.com/article/10.1007/s00023-023-01325-x

#bibliography(
  title: "References",
  style: "ieee",
  "bibliography.bib"
)
