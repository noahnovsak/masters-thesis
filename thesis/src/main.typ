#import "conf.typ": conf
#import "@preview/drafting:0.2.2": *
#import "@preview/algorithmic:1.0.7"
#import "@preview/cheq:0.3.0": checklist
#import algorithmic: *

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

#let PPT2 = $"PPT"^2$

#let make-statement(kind, label: none) = {
  let display-label = if label == none {
    kind
  } else {
    label
  }

  (title: none, body) => context {
    let chapter = query(heading.where(level: 1, outlined: true).before(here())).len()

    // Use a chapter-scoped counter so numbering keeps working even if heading numbers are hidden.
    let c = counter(kind + "-" + str(chapter))
    c.step()

    block(
      above: 1em,
      below: 1em,
      [
        #strong([#display-label #chapter.#c.display("1.1")])#if title != none {
          [ (#emph[#title]).]
        } else {
          [.]
        }
        #h(0.5em)
        #body
      ],
    )
  }
}

#let conjecture = make-statement("conjecture", label: "Conjecture")
#let theorem = make-statement("theorem", label: "Theorem")
#let proposition = make-statement("proposition", label: "Proposition")
#let corollary = make-statement("corollary", label: "Corollary")
#let definition = make-statement("definition", label: "Definition")
#let remark = make-statement("remark", label: "Remark")
#let definition = make-statement("definition", label: "Definition")

#let reference = margin-note("reference")

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

We shall work with bipartite states on the space $cal(H)_A times.o cal(H)_B$ (with dimensions $m, n$).

All PPT maps on $M_2(CC)$ are separable due to the Peres-Horodecki separation criterion @Peres_1996.

The PPT2 conjecture holds trivially for $n = 2$, and has been proven for $n eq 3$ @Chen_2019. Additionally it has been proven to hold in all dimensions for maps with specific (physically relevant) properties @Singh_2022. 

Obviously the composition of two PPT maps is still PPT.

A quntum state is a vector $lr(| Psi chevron.r)$ in $CC^n$, also represented as a positive semidefinite (PSD) density matrix $rho = lr(| Psi chevron.r chevron Psi |)$. A quantum channel is a completely positive trace-preserving (CPTP) linear map $Phi: cal(M)_n -> cal(M)_m$.
#inline-note[Pure state/Mixed state/ density matrix definition needed?]

The set of all quantum states is exactly the cone of positive semidefinite matrices, and the set of all quantum channels is exactly the cone of completely positive trace-preserving maps. The Choi-Jamiolkowski isomorphism provides a convenient way to represent linear maps as matrices, which allows us to apply semidefinite programming techniques for testing properties like complete positivity and entanglement.

A state $rho$ in $cal(H)_A times.o cal(H)_B$ is _seperable_ if it can be written as a convex combination of product states, i.e. $rho = sum_i p_i rho_i^A times.o rho_i^B$ with $p_i >= 0$, $sum_i p_i = 1$, and $rho_i^A in cal(H)_A$, $rho_i^B in cal(H)_B$. Otherwise it is _entangled_.

Determining whether a given state is separable or entangled is a hard problem in general @Gharibian_2009.
#inline-note[A necessary condition for separability is that the state has a positive partial transpose (PPT), but this condition is not sufficient in higher dimensions.]
As stated previously, a separable state can be written as a linear combination of product states, i.e. $rho = rho_A times.o rho_B$. Since $rho_A$ and $rho_B$ are both positive semidefinite matrices the partial transpose $rho^Gamma = rho_A times.o rho_B^T$ is also positive semidefinite. This means that all separable states are necessarily PPT. However, the converse is not always true; there are entangled states that are also PPT when $m n >= 6$@Horodecki_2009.

A linear map $Phi: M_n -> M_m$ is: _positive_ if it maps positive semidefinite matrices to positive semidefinite matrices, _k-positive_ if $I_k times.o Phi$ is a positive map, and _completely positive_ if it is $k$-positive for all $k$.

Let $M^+$ denote the cone of positive semidefinite matrices.

#definition(title: "Positive map")[
  A linear map $Phi: M_n -> M_m$ is positive if $X in M_n^+ => Phi(X) in M_m^+$.
]

#definition(title: "k-positive map")[
  A linear map $Phi: M_n -> M_m$ is k-positive if $I_k times.o Phi$ is positive.
]

#definition(title: "Completely positive map")[
  A linear map $Phi: M_n -> M_m$ is completely positive if it is $k$-positive for all $k$.
]

#definition(title: "Quantum state")[
  A quantum state is a vector $| Psi >$ in $CC^n$.
]

#definition(title: "Quantum channel")[
  A quantum channel is a completely positive trace-preserving (CPTP) linear map $Phi: cal(M)_n -> cal(M)_m$.
]

#definition(title: "Choi-Jamiolkowski isomorphism")[
  Given a linear map $Phi: cal(M)_n -> cal(M)_m$, we call $ C_Phi = (I times.o Phi)(w) = sum_(i j) E_(i j) times.o Phi(E_(i j)) $ the Choi matrix of $Phi$. The map $Phi arrow.bar C_Phi$ defines an isomorphism called the Choi-Jamiolkowski isomorphism between linear maps $Phi: cal(M)_n -> cal(M)_m$ and matrices in $cal(M)_n times.o cal(M)_m$.
]
Under the Choi-Jamiolkowski isomorphism completely positive maps correspond to positive semidefinite matrices.

#definition(title: "Partial transpose")[
  Given a state $rho = sum rho_A times.o rho_B$ acting on $cal(H)_A times.o cal(H)_B$, its partial transpose (with respect to subsystem $B$) is defined as $rho^(T_B) = (I times.o T)(rho) = sum rho_A times.o rho_B^T$.
]
In this context we are almost always working with symmetric matrices, where the choice of subsystem for partial transpose is irrelevant, so we simply write $rho^Gamma$ for the partial transpose.

A state $rho$ is said to be PPT (positive partial transpose) if its partial transpose $rho^Gamma$ is positive semidefinite.

#inline-note[In this work, we move back and forth between discussingstates, channels, matrices, and maps. Due to the Choi-Jamiolkowski isomorphism, these are all essentially the same objects in different representations. Their properties (like positivity, entanglement, etc.) can be checked in any representation, but some operations (like composition) are more naturally expressed in the map representation, while others (like partial transpose) are more naturally expressed in the matrix representation. When working with the numerical implementation, we always use the matrix representation, but we keep the map representation in mind for conceptual clarity and to ensure that our operations correspond to the intended mathematical objects.]

A quantum channel is a completely positive trace-preserving (CPTP) linear map $Phi: cal(M)_n -> cal(M)_m$. We say it is entanglement breaking/PPT/etc. if its Choi matrix has the corresponding property. 

#definition[
  PPT (positive partial transpose)
]

#definition[
  Entanglement
]

Any positive but not completely positive map $Phi: cal(H)_A -> cal(H)_B$ provides the necessary separability criterion $(I times.o Phi)(rho) succ.eq 0$ @Horodecki_2009 @Horodecki_1996. In the case of the transposition map $T$ we get the well-known PPT criterion $(I times.o T)(rho) succ.eq 0$ @Peres_1996 @Horodecki_1996.

#definition[
  Entanglement witnesses are observables that completely characterize the set of separable states and allow us to detect entanglement physically. For every entangled state $rho$, there exists a Hermitian operator $W$ such that $tr(W rho) < 0$ and $tr(W sigma) >= 0$ for all separable states $sigma$. Such an operator $W$ is called an entanglement witness for the state $rho$.
]

The positive but not completely positive (PNCP) maps and entanglement witnesses are linked by the Choi-Jamiolkowski isomorphism $C_Phi = (I times.o Phi)(omega)$. However, an important observation is that while $tr(C_Phi rho) >= 0$ is equivalent to $(I times.o Phi)(rho) succ.eq 0$, a particular witness is not equivalent to a positive map associated via isomorphism: the map proves a stronger condition @Horodecki_2009.

Decomposable witnesses can be written as $W = P + Q^Gamma$, where $P$, $Q$ are positive operators @Lewenstein_2000. Such witnesses cannot detect PPT-entangled states, so we are particularly interested in non-decomposable witnesses.

=== Semidefinite programming and DPS hierarchy
Consider a quantum state $rho$ on systems $cal(H)_A$ and $cal(H)_B$, with dimensions $m$ and $n$. It is separable if there exist pure states $lr({ | phi_i chevron.r_A })$, $lr({ | chi_i chevron.r_B })$ and probabilities ${p_i}$ such that $rho = sum_i p_i lr(|phi_i chevron.r chevron phi_i| times.o |chi_i chevron.r chevron chi_i |)$.

Let $S E P_(n,k) := "conv"lr({|phi_1 chevron.r chevron phi_1| times.o ... times.o |phi_k chevron.r chevron phi_k| : |phi_1 chevron.r, ..., |phi_k chevron.r in B(CC^n)})$ denote the set of separable states, where conv is the convex hull and _B_ is the set of unit vectors.

The set of separable states can be described as $S E P_n = lr("cone"{x x^* times.o y y^* : x, y in CC^n, ||x|| = ||y|| = 1} subset.eq H(CC^n times.o CC^n))$.

Testing entanglement can be reduced to the problem $h_(S E P)(M) = lr(max{"Tr"(M rho) : rho in S E P})$.

Consider a separable state $rho in H(CC^n times.o CC^n) := sum_i lambda_i x_i x_i^* times.o y_i y_i^*$. Then, for integers $k >= 1$ the extended states $rho_k = sum_i lambda_i x_i x_i^* times.o (y_i y_i^*)^(times.o k)$, have the following properties: $ (I_n times.o Pi_k) rho_k (I_n times.o Pi_k) = rho_k, \ rho_k^Gamma_s succ.eq 0 space forall space s = 0, 1, ..., k, \ "Tr"_([2:k])(rho_k) = rho. $

These states form a hierarchy of semidefinite relaxations that approximate the cone of separable states $S E P_n$ @Doherty_2004 $ "DPS"_n^k = {rho in H(CC^n times.o CC^n) : exists space rho_k in H(CC^n times.o (CC^n)^(times.o k))}. $

Here we note  few properties of the hierarchy; 1. $"DPS"_n^1$ is equivalent to the PPT criterion, 2. $S E P_n subset.eq "DPS"_n^k$ and $"DPS"_n^(k+1) subset.eq "DPS"_n^k$ for all $k >= 1$, 3. the hierarchy is complete $inter.big_(k >= 1) "DPS"_n^k = S E P_n$ @Doherty_2004.

A semidefninite program (SDP) is a type of convex optimization problem typically expressed as $ "minimize"& c^T bold(x) \ "subject to"& F(bold(x))  succ.eq 0, $ where $c$ is a given vector, $F(bold(x)) = F_0 + sum_i x_i F_i$ for some fixed hermitian matrices $F_i$, and $bold(x) = (x_1, ..., x_n)$ is the vector over which the optimization is performed. Importantly, for any such SDP (called the _primal_ problem) there exists an associated _dual_ problem $ "maximize"& -"Tr"[F_0 Z] \ "subject to"& Z  succ.eq 0 \ &"Tr"[F_i Z] = c_i, $ where $Z$ is the hermitian optimization matrix. When $c = 0$ we call the primal problem a _feasibility_ problem, and the dual problem a _certificate_ problem. In this case, if the primal is infeasible, the dual is feasible and provides a certificate of infeasibility.

In the case of the $"DPS"_n^t$ hierarchy, whenever the primal problem is infeasible (a PPT symmetric extension of $rho$ cannot be found), the certificate provided by the dual problem is an entanglement witness for the state $rho$.

Much work has been done recently on semidefinite hierarchies for states with specific properties, i.e. diagonal unitary invariant bipartite
quantum states @Britz_2025, the trouble here is that we cannot use them for our purposes since the PPT2 conjecture has already been proven for this family of states @Singh_2022, so we need to work with the more general DPS hierarchy.

However, even in the general case some improvements have been made, for example adding KKT constraints to the SDP to achieve finite convergence @Harrow_2017. The downside here is that an entanglement witness cannot be as cleanly extracted from the dual solution, so verification becomes more complex.

== PPT2 Conjecture
#conjecture[
  The composition of two PPT linear maps is entanglement breaking.
]

Let $M_n$ denote complex $n x n$ matrices and $M_n^+$ the cone of positive semidefinite matrices.

#definition(title: "Positive map")[
  A linear map $Phi: M_n -> M_m$ is positive if $X in M_n^+ => Phi(X) in M_m^+$.
]

#definition(title: "k-positive and completely positive map")[
  For $k >= 1$, map $Phi$ is $k$-positive if $I_k times.o Phi$ is positive. It is completely positive (CP) if it is $k$-positive for all $k$ @Chen_2019 @jin2020.
]

#inline-note[Although channels are trace preserving by definition, the computational tests in this thesis are primarily positivity/sign based. Therefore, intermediate computations may use non-normalized maps, then normalize where needed for interpretation.]

#inline-note[a quantum state is usually represented by a positive semidefinite operator normalized with unit trace. In this work, we mainly work with unnormalized quantum states and consider its convex cone.]

#inline-note[A well-known necessary condition for a state $rho_(A B)$ to be separable is that is has a positive partial transpose (PPT).]

#definition(title: "PPT map")[
  A linear map $Phi: M_n -> M_m$ is PPT if $(I_n o T_m)(C_Phi) >= 0$, where $T_m$ is the transpose map on $M_m$.
]

#definition(title: "Entanglement-breaking map")[
  A map $Phi: M_n -> M_n$ is entanglement breaking if $(I_n compose Phi)(rho)$ is separable for every bipartite state $rho$.
]

#conjecture(title: "PPT2 conjecture")[
  If $Phi$ is a PPT map, then $Phi compose Phi$ is entanglement breaking @Christandl_2019. An equivalent formulation is that the composition of any two PPT maps is entanglement breaking.
]

== Entanglement and Separability
#definition(title: "Bipartite state")[
  A bipartite state is $rho in M_n times.o M_m$ with $rho >= 0$ and $"tr"(rho)=1$.
]

#definition(title: "Separable state")[
  State $rho$ is separable if $rho = sum_i p_i rho_i^A x rho_i^B$ with $p_i >= 0$, $sum_i p_i = 1$, and local states $rho_i^A, rho_i^B$. Otherwise it is entangled @Horodecki_2009.
]

#definition(title: "Partial transpose and PPT state")[
  Let $Gamma_B := id_n o T_m$. A state is PPT if $Gamma_B(rho) >= 0$ @Peres_1996 @Horodecki_1996.
]

#proposition(title: "PPT criterion, necessity")[
  Every separable state is PPT @Peres_1996 @Horodecki_1996.
]

#remark[
  In low dimensions (notably $2 x 2$ and $2 x 3$), PPT is also sufficient for separability. In higher dimensions this fails, and PPT-entangled states exist @Horodecki_2009.
]

#inline-note[Witnesses are one-sided certificates: a valid negative expectation certifies entanglement, while non-negativity does not certify separability.]

#definition(title: "DPS hierarchy, informal")[
  For each level $ell$, DPS defines an SDP relaxation set $S_ell$ based on symmetric extensions. If $rho in.not S_ell$, entanglement is certified. As $ell$ increases, relaxations tighten @Doherty_2004.
]

The hierarchy is complete in the limit and is a standard baseline in practice. Improved variants exist, but finite-level trade-offs remain important @Harrow_2017.

== Positive Maps and Choi Representation
#definition(title: "Choi matrix")[
  For $Phi: M_n -> M_n$, define
  $
    C_Phi = (I_n times.o Phi)(Omega Omega^*),
  $
  where $Omega = sum_(i=1)^n e_i x e_i$.
]

#theorem(title: "Choi criterion")[
  Map $Phi$ is completely positive iff $C_Phi >= 0$ @Choi_1975.
]

#definition(title: "PPT in Choi form")[
  Map $Phi$ is PPT iff $C_Phi >= 0$ and $C_Phi^Gamma >= 0$.
]

For composition, Choi matrices are combined via ampliation rather than ordinary matrix multiplication:
$
  C_(Phi compose Psi) = (I_n times.o C_Phi)(C_Psi).
$

== PNCP Maps
Positive but not completely positive (PNCP) maps are useful because they induce witness families through Choi duality. This thesis uses PNCP construction as a complementary route to DPS-based witness extraction, with SOS checks and post-solver revalidation @Fang_2020 @bhardwaj2020 @phdthesis.

Given two matrix spaces $M_n$ and $M_m$, a linear map $Phi: M_n -> M_m$ is _positive_ if $Phi(A) succ.eq 0 forall A succ.eq 0$.
For a given $k in NN$ such maps induce the _ampliation_ $Phi^k = I_k times.o Phi: M_k times.o M_n -> M_k times.o M_m$. If $Phi^k$ is positive, we say that $Phi$ is _$k$-positive_. If $Phi$ is $k$-positive for all $k$, we say that it is _completely positive_ (CP) @bhardwaj2020.

Further, each linear map $Phi: M_n -> M_m$ corresponds to a quadratic polynomial $p_Phi (bold(x),bold(y)) in RR[bold(x), bold(y)]_(n,m) := bold(y)^T Phi(bold(x) bold(x)^T) bold(y)$. It is known that $Phi$ is positive when $p_Phi$ is non-negative, and completely positive when $p_Phi$ is a sum of squares (SOS) on $RR^(n+m)$ @Klep_2017.

#algorithm-figure(
  "PNCP construction",
  {
    Line[Generate random points $x in RR^n$, $y in RR^m$.]
    Line[Use $x$, $y$ to create bilinear forms ${h_0, ..., h_t}$]
    Line[Generate $f in.not lr(chevron h_0, ..., h_t chevron.r)$ such that $f != text(S O S)$.]
    Line[Find $delta$ small enough that $F_delta = delta f + h_0^2 + ... + h_t^2 >= 0$.]
  },
) <pncp-algo>

The first three steps are purely linear algebra, and the last step is a semidefinite program. The resulting polynomial $F_delta$ is non-negative and non-SOS, so the corresponding map is positive but not completely positive @Klep_2017.

Our constructions are intrinsically numerical, and the resulting maps are represented as floating-point matrices. This is sufficient for quick tests and exploration, but numerical noise can lead to false positives near the boundary of positivity. Therefore, we apply a post-solver rationalization and re-validation step to confirm any positive detections.

The semidefinite programs arising from our SOS relaxations are feasibility problems of the form $ G succ.eq 0 \ s.t. space chevron A_i, G chevron.r = b_i, quad i = 1, ..., m $ where $A_i$, $b_i$ are obtained from the problem data.

Let $G$ be a positive definite feasible point satisfying $mu = min("eig"(G)) > ||(chevron A_i, G chevron.r) - b_i|| = epsilon$, then a rational feasible point $hat(G)$ can be obtained by computing a rational approximation $tilde(G)$ with $|| G - tilde(G) ||^2 + epsilon^2 <= mu^2$, then projecting it back to the affine space defined by the constraints $chevron A_i, G chevron.r = b_i$. 

Steps 1-3 can be performed over $QQ$, leading to rational forms $h_j$, $f$. However in step 4 we need to solve an SDP to find a suitable $delta$, which will be floating-point. So instead of rationalizing the entire problem beforehand, we only rationalize the solution after solving the SDP. From our construction we know that each feasible $G$ will have at least an $e$-dimensional nullspace @Klep_2017. Therefore, we can create a rationalized $tilde(G) = mat(tilde(G)_11, tilde(G)_12; tilde(G)_12^T, tilde(G)_22)$, where we set $tilde(G)_11$ and $tilde(G)_12$ both equal to $0$, then rationalize the remaining coefficients $tilde(G)_22$. Finally we run an additional feasibility SDP to verify $tilde(G) succ.eq 0$.

= Methods
This chapter defines the computational workflow used in the current implementation.

== Chapter Goals and Design Constraints
The goal is reproducibility with mathematical traceability: each code step corresponds to a defined object or operation from Chapter 2.

Design constraints:
1. separability testing difficulty in relevant dimensions,
2. numerical fragility near SDP feasibility boundaries,
3. need to separate broad exploration from strict validation.

== Inputs and Outputs
Inputs are dimensions $(n,m)$, random seeds, and solver parameters. Outputs are candidate states/maps, composed Choi objects, witness scores, and validation artifacts.

The workflow stores intermediate objects to enable replay, debugging, and independent certificate checks.

== Pipeline A: Randomized PPT Composition + DPS
#algorithm-figure(
  "Randomized PPT Composition + DPS",
  {
    Procedure(
      "Pipeline-A", ("D", "N_d", "profile_DPS"),
      {
        Line[*Input:* dimension set $D$; sample counts $N_d$; DPS level.]
        Line[*Output:* screened candidate list with DPS-derived evidence.]
        LineBreak
        For(
          $d in D$,
          {
            Line[Generate $N_d$ random candidates under a rank policy.]
            Line[Filter by PSD and PPT constraints.]
            Line[Convert to common Choi representation.]
            Line[Build self and pairwise compositions via ampliation.]
            Line[Run DPS tests and witness evaluations.]
            Line[Persist solver logs and intermediate artifacts.]
          },
        )
      },
    )
  },
)<algo-1>

Pipeline A is high-throughput and broad-coverage; positive detections are treated as preliminary until stricter validation.

== Pipeline B: Polynomial/SOS Witness Construction
#algorithm-figure(
  "Polynomial/SOS Witness Construction",
  {
    Procedure(
      "Pipeline-B",
      ("params_map", "level_SOS", "tau_eig", "targets_A"),
      {
        Line[*Input:* map parameters; SOS level; eigenvalue policy; targets.]
        Line[*Output:* validated witness candidates and evaluation logs.]
        LineBreak
        Line[Construct PNCP candidates via polynomial parametrization.]
        Line[Solve SOS programs and extract Gram matrices.]
        Line[Clean near-zero spectrum and reconstruct certificates.]
        Line[Re-test cleaned candidates in a fresh feasibility pass.]
        Line[Apply surviving witnesses to composed targets.]
      },
    )
  },
)<algo-2>

Pipeline B is lower-throughput and higher-scrutiny, aimed at reducing false positives from floating-point noise.

== Composition Operator and Ampliation
The implementation in `code/src/ppt2.jl` computes composition with explicit block conversions (`mat2block`, `block2mat`) and the `ampliation` operator. This realizes the Choi composition rule directly and keeps basis handling explicit.

== Rationalization and Certificate Validation
For SOS outputs, the implementation applies a post-solver validation loop:
1. extract Gram representation,
2. clean near-zero eigencomponents,
3. reconstruct polynomial candidate,
4. rerun feasibility on the reconstructed object.

This is implemented in `solve_sos` through the `zero_g` branch and follow-up checks.

== Numerical Reliability Protocol
Evidence is accepted only after layered checks:
1. initial detection (DPS or witness),
2. independent cross-check where available,
3. post-processed certificate re-validation.

Tiny negative values are treated as inconclusive unless they remain stable under the full validation chain.

== Implementation Architecture
The active implementation is in the `code` directory, centered on `code/src/ppt2.jl`.

Exported primitives:
1. `rand_ppt(n,m)` for candidate generation,
2. `ampliation(A,B,n,m)` for composition,
3. `pncp_mat(n,m)` for PNCP/SOS witness route.

Supporting internals include `kernel_basis`, `quadratic_form`, `pncp_algorithm`, `solve_sos`, and representation transforms (`poly2mat`, `mat2block`, `block2mat`).

== Complexity and Practical Limits
Both routes rely on SDPs whose size grows quickly with dimension and relaxation depth. Practical bottlenecks are:
1. memory growth,
2. solver instability near degeneracy,
3. long-tail runtime on difficult instances.

Therefore, the architecture separates broad search from strict validation and preserves artifacts for selective reruns.

= Experimental Evaluation
This manuscript version emphasizes mathematical and methodological development. Experimental outputs are treated as supporting context for the workflow design rather than the main argumentative core.

= Methodological Limitations
The approach remains exploratory and has clear limits:
1. finite DPS levels may miss hard PPT-entangled states,
2. generated witness families are not guaranteed to be optimal,
3. boundary cases remain numerically delicate,
4. high dimensions and hierarchy levels can become computationally prohibitive.

= Conclusion
The thesis establishes a reproducible computational research program for the #PPT2 conjecture by aligning formal map-theoretic formulations with implementation-level operators and validation procedures.

Its principal value at this stage is methodological: it defines a transparent path from candidate construction to certificate verification, enabling future counterexample searches and stronger structural analyses.

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
