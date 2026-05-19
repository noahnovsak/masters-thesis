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
This chapter formalizes notation and statements used by the implementation.

== Preliminaries
inner product: $"Tr"[A B] = chevron A, B chevron.r$

Let $M_n$ denote the set of complex $n times n$ matrices and $M_n^+$ the cone of positive semidefinite matrices.

We shall work with bipartite states on the space $cal(H)_A times.o cal(H)_B$ (with dimensions $m, n$).

All PPT maps on $M_2(CC)$ are separable due to the Peres-Horodecki separation criterion @Peres_1996. The #PPT2 conjecture holds trivially for $n = 2$, and has been proven for $n eq 3$ @Chen_2019. Additionally it has been proven to hold in all dimensions for maps with specific (physically relevant) properties @Singh_2022. The composition of two PPT maps is trivially PPT.

The set of all quantum states is exactly the cone of positive semidefinite matrices, and the set of all quantum channels is exactly the cone of completely positive trace-preserving maps. The Choi-Jamiolkowski isomorphism provides a convenient way to represent linear maps as matrices, which allows us to apply semidefinite programming techniques for testing properties like complete positivity and entanglement.

A state $rho$ in $cal(H)_A times.o cal(H)_B$ is _separable_ if it can be written as a convex combination of product states, i.e. $rho = sum_i p_i rho_i^A times.o rho_i^B$ with $p_i >= 0$, $sum_i p_i = 1$, and $rho_i^A in cal(H)_A$, $rho_i^B in cal(H)_B$. Otherwise it is _entangled_.

#definition(name: "Entanglement-breaking map")[
  A map $Phi: M_n -> M_n$ is entanglement breaking if $Phi(rho)$ is separable for every bipartite state $rho$.
]

Determining whether a given state is separable or entangled is a hard problem in general @Gharibian_2009. As stated, a separable state can be written as a linear combination of product states, i.e. $rho = rho_A times.o rho_B$. Since $rho_A$ and $rho_B$ are both positive semidefinite matrices the partial transpose $rho^Gamma = rho_A times.o rho_B^T$ is also positive semidefinite. This means that all separable states are necessarily PPT. However, the converse is not always true; there are entangled states that are also PPT when $m n > 6$ @Horodecki_2009.

A linear map $Phi: M_n -> M_m$ is: _positive_ if it maps positive semidefinite matrices to positive semidefinite matrices, _k-positive_ if $I_k times.o Phi$ is a positive map, and _completely positive_ if it is $k$-positive for all $k$.

#definition(name: "Positive map")[
  A linear map $Phi: M_n -> M_m$ is positive if $X in M_n^+ => Phi(X) in M_m^+$.
]

#definition(name: "k-positive and completely positive map")[
  For $k >= 1$, map $Phi$ is $k$-positive if $I_k times.o Phi$ is positive. It is completely positive (CP) if it is $k$-positive for all $k$ @Chen_2019 @jin2020.
]

Given two matrix spaces $M_n$ and $M_m$, a linear map $Phi: M_n -> M_m$ is _positive_ if $Phi(A) succ.eq 0 space forall space A succ.eq 0$.
For a given $k in NN$ such maps induce the _ampliation_ $Phi^k = I_k times.o Phi: M_k times.o M_n -> M_k times.o M_m$. If $Phi^k$ is positive, we say that $Phi$ is _$k$-positive_. If $Phi$ is $k$-positive for all $k$, we say that it is _completely positive_ (CP) @bhardwaj2020.

#definition(name: "Choi-Jamiolkowski isomorphism")[
  Given a linear map $Phi: M_n -> M_m$, we call $ C_Phi = (I times.o Phi)(w) = sum_(i j) E_(i j) times.o Phi(E_(i j)) $ the Choi matrix of $Phi$. The map $Phi arrow.bar C_Phi$ defines an isomorphism called the Choi-Jamiolkowski isomorphism between linear maps $Phi: M_n -> M_m$ and matrices in $M_n times.o M_m$.
]
Under the Choi-Jamiolkowski isomorphism completely positive maps correspond to positive semidefinite matrices (PSD) @Choi_1975.

#definition(name: "Partial transpose")[
  Given a state $rho = sum rho_A times.o rho_B$ acting on $cal(H)_A times.o cal(H)_B$, its partial transpose (with respect to subsystem $B$) is defined as $rho^(T_B) = (I times.o T)(rho) = sum rho_A times.o rho_B^T$.
]
In this context we are almost always working with symmetric matrices, where the choice of subsystem for partial transpose is irrelevant, so we simply write $rho^Gamma$ for the partial transpose.

A state $rho$ is said to be PPT (positive partial transpose) if its partial transpose $rho^Gamma$ is positive semidefinite.

#definition(name: "PPT map")[
  A linear map $Phi: M_n -> M_m$ is PPT if $(I_n compose T_m)(C_Phi) succ.eq 0$, where $T_m$ is the transpose map on the second subsystem.
]

#definition(name: "PPT in Choi form")[
  Map $Phi: M_n -> M_m$ is PPT iff $C_Phi succ.eq 0$ and $C_Phi^Gamma succ.eq 0$.
]

A quantum channel is a completely positive trace-preserving (CPTP) linear map $Phi: M_n -> M_m$. Due to the Choi-Jamiolkowski isomorphism, states, channels, matrices, and maps are essentially equivalent representations of the same objects. Their properties (like positivity and entanglement) can be checked in any representation, but some operations (like composition) are more naturally expressed in the map representation, while others (like partial transpose) are more naturally expressed in the matrix representation. 

#definition(name: "Entanglement witness")[
  Entanglement witnesses are observables that completely characterize the set of separable states and allow us to detect entanglement physically. For every entangled state $rho$, there exists a Hermitian operator $W$ such that $tr(W rho) < 0$ and $tr(W sigma) >= 0$ for all separable states $sigma$. Such an operator $W$ is called an entanglement witness for the state $rho$.
]

Positive but not completely positive (PNCP) maps and entanglement witnesses are linked by the Choi-Jamiolkowski isomorphism $C_Phi = (I times.o Phi)(omega)$. However, an important observation is that while $tr(C_Phi rho) >= 0$ is equivalent to $(I times.o Phi)(rho) succ.eq 0$, a particular witness is not equivalent to a positive map associated via isomorphism: the map proves a stronger condition @Horodecki_2009.

Decomposable witnesses can be written as $W = P + Q^Gamma$, where $P$, $Q$ are positive operators @Lewenstein_2000. Such witnesses cannot detect PPT-entangled states, so we are particularly interested in non-decomposable witnesses.

#inline-note[Witnesses are one-sided certificates: a valid negative expectation certifies entanglement, while non-negativity does not certify separability.]

When computing the composition of two maps numerically, we are in reality working with their Choi matrices. Therefore, we cannot simply multiply them, instead we compute:
$
  C_(Phi compose Psi)
  = sum_(i,j) E_(i j) times.o Phi(C^Psi_(i j))
  = sum_(i,j) E_(i j) times.o sum_(k,l) (C^Psi_(i j))_(k l) Phi(E_(k l))
  = sum_(i,j,k,l,p,q) (C^Psi_(i j))_(k l) (C^Phi_(k l))_(p q) E_(i j) times.o E_(p q)
$
$ C^(Phi compose Psi)_(i p, j q) = sum_(k,l) C^Psi_(i k,p l) C^Phi_(k j, l q). $

== Semidefinite programming and DPS hierarchy
Consider a quantum state $rho$ on systems $cal(H)_A$ and $cal(H)_B$, with dimensions $m$ and $n$. It is separable if there exist pure states $lr({ | phi_i chevron.r_A })$, $lr({ | chi_i chevron.r_B })$ and probabilities ${p_i}$ such that $rho = sum_i p_i lr(|phi_i chevron.r chevron phi_i| times.o |chi_i chevron.r chevron chi_i |)$.

Let us denote the set of all separable states $S E P_(n,k) := "conv"lr({|phi_1 chevron.r chevron phi_1| times.o ... times.o |phi_k chevron.r chevron phi_k| : |phi_1 chevron.r, ..., |phi_k chevron.r in B(CC^n)})$, where conv is the convex hull and _B_ is the set of unit vectors.

Testing entanglement can be reduced to the problem $h_(S E P)(M) = lr(max{"Tr"(M rho) : rho in S E P})$ @Harrow_2017. #margin-note[expand]

Consider a separable state $rho in H(CC^n times.o CC^n) := sum_i lambda_i x_i x_i^* times.o y_i y_i^*$. Then, for integers $k >= 1$ the extended states $rho_k = sum_i lambda_i x_i x_i^* times.o (y_i y_i^*)^(times.o k)$, have the following properties: $ (I_n times.o Pi_k) rho_k (I_n times.o Pi_k) = rho_k, \ rho_k^Gamma_s succ.eq 0 space forall space s = 0, 1, ..., k, \ "Tr"_([2:k])(rho_k) = rho. $

These states form a hierarchy of semidefinite relaxations that approximate the cone of separable states $S E P_n$ @Doherty_2004 $ "DPS"_n^k = {rho in H(CC^n times.o CC^n) : exists space rho_k in H(CC^n times.o (CC^n)^(times.o k))}. $

Here we note  few properties of the hierarchy; 1. $"DPS"_n^1$ is equivalent to the PPT criterion, 2. $S E P_n subset.eq "DPS"_n^k$ and $"DPS"_n^(k+1) subset.eq "DPS"_n^k$ for all $k >= 1$, 3. the hierarchy is complete $inter.big_(k >= 1) "DPS"_n^k = S E P_n$ @Doherty_2004.

A semidefninite program (SDP) is a type of convex optimization problem typically expressed as $ "minimize"& c^T bold(x) \ "subject to"& F(bold(x))  succ.eq 0, $ where $c$ is a given vector, $F(bold(x)) = F_0 + sum_i x_i F_i$ for some fixed hermitian matrices $F_i$, and $bold(x) = (x_1, ..., x_n)$ is the vector over which the optimization is performed. Importantly, for any such SDP (called the _primal_ problem) there exists an associated _dual_ problem $ "maximize"& -"Tr"[F_0 Z] \ "subject to"& Z  succ.eq 0 \ &"Tr"[F_i Z] = c_i, $ where $Z$ is the hermitian optimization matrix. When $c = 0$ we call the primal problem a _feasibility_ problem, and the dual problem a _certificate_ problem. In this case, if the primal is infeasible, the dual is feasible and provides a certificate of infeasibility.

In the case of the $"DPS"_n^t$ hierarchy, whenever the primal problem is infeasible (a PPT symmetric extension of $rho$ cannot be found), the certificate provided by the dual problem is an entanglement witness for the state $rho$.

Much work has been done recently on semidefinite hierarchies for states with specific properties, i.e. diagonal unitary invariant bipartite quantum states @Britz_2025, the trouble here is that we cannot use them for our purposes since the PPT2 conjecture has already been proven for this family of states @Singh_2022, so we need to work with the more general DPS hierarchy.

However, even in the general case some improvements have been made, for example adding KKT constraints to the SDP to achieve finite convergence @Harrow_2017. The downside here is that an entanglement witness cannot be as cleanly extracted from the dual solution, so verification becomes more complex.

== PNCP Maps and Sum-of-Squares Polynomials
Positive but not completely positive (PNCP) maps are useful for two reasons: every entanglement witness corresponds to a PNCP map, and they can be used to construct bound entangled states, i.e. entangled states with positive partial transpose (PPT). In this thesis we use PNCP maps as an alternative route to DPS-based approaches for testing entanglement and witness extraction, enabling exact artithmetic validation @Fang_2020 @bhardwaj2020 @phdthesis.

Each linear map $Phi: M_n -> M_m$ corresponds to a quadratic polynomial $p_Phi (bold(x),bold(y)) in RR[bold(x), bold(y)]_(n,m) := bold(y)^T Phi(bold(x) bold(x)^T) bold(y)$. It is known that $Phi$ is positive if and only if $p_Phi$ is non-negative, and completely positive if and only if $p_Phi$ is a sum of squares (SOS) on $RR^(n+m)$ @Klep_2017.

So now, in addition to matrices and maps we have a third, polynomial, representation of the same objects. This allows us to use tools from polynomial optimization, in particular SOS relaxations. @kmsz provides a way to construct PNCP maps by generating random polynomials that are non-negative but not SOS.

#algorithm-figure(
  "KMSZ construction for PNCP maps",
  {
    Line[Generate random points $x in RR^n$, $y in RR^m$.]
    Line[Use $x$, $y$ to create bilinear forms ${h_0, ..., h_t}$]
    Line[Generate $f in.not lr(chevron h_0, ..., h_t chevron.r)$ such that $f != text(S O S)$.]
    Line[Find $delta$ small enough that $F_delta = delta f + h_0^2 + ... + h_t^2 >= 0$.]
  },
) <kmsz>

The first three steps are purely linear algebra, and the last step is a semidefinite program, which we will describe in more detail later. The resulting polynomial $F_delta$ is non-negative and non-SOS, so the corresponding map is positive but not completely positive @Klep_2017.

=== Rationalization
Steps 1-3 of @kmsz can be performed over $QQ$, leading to rational forms $h_j$, $f$. However in step 4 we need to solve an SDP to find a suitable $delta$. Any state-of-the-art solver will produce a floating-point solution. So instead of rationalizing the entire problem beforehand, we only rationalize the solution after solving the SDP.

The semidefinite programs arising from our SOS relaxations are feasibility problems of the form $ G succ.eq 0 \ s.t. space chevron A_i, G chevron.r = b_i, quad i = 1, ..., m $ where $A_i$, $b_i$ are obtained from the problem data.

Let $G$ be a positive definite feasible point satisfying $mu = min("eig"(G)) > ||(chevron A_i, G chevron.r) - b_i|| = epsilon$, then a rational feasible point $hat(G)$ can be obtained by computing a rational approximation $tilde(G)$ with $|| G - tilde(G) ||^2 + epsilon^2 <= mu^2$, then projecting it back to the affine space defined by the constraints $chevron A_i, G chevron.r = b_i$. 

From our construction we know that each feasible $G$ will have at least an $e$-dimensional nullspace @Klep_2017. Therefore, we can create a rationalized $tilde(G) = mat(tilde(G)_11, tilde(G)_12; tilde(G)_12^T, tilde(G)_22)$, where we set $tilde(G)_11$ and $tilde(G)_12$ are both equal to $0$, then rationalize the remaining coefficients $tilde(G)_22$. Finally we run an additional feasibility SDP to verify $tilde(G) succ.eq 0$.

== PPT2 Conjecture
#theorem(name: "PPT2 Conjecture")[
  Given a PPT map $phi$, the composition $phi compose phi$ is separable.@Christandl_2019
]<ppt2>

Equivalently to @ppt2, if $phi_1$, $phi_2$ are two PPT maps, then $phi_1 compose phi_2$ is separable @Chen_2019.

= Methods
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
