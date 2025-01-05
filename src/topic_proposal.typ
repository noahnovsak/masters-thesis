#set document(
    title: [Masters Thesis Topic Proposal],
    author: "Noah Novšak",
    date: datetime.today(),
)
#set par(
    justify: true,
    leading: 0.6em,
    spacing: 1.8em,
)
#set page(
    numbering: "1",
    margin: 2.5cm,
    footer-descent: 50%,
)
#set text(
    font: "New Computer Modern",
)

#let PPT2 = $"PPT"^2$


Noah Novšak \
Povšetova ulica 44, 1000 Ljubljana \
Computer and Information Science, MAG \
63230470

*Committee for Student Affairs* \
University of Ljubljana, Faculty of Computer and Information Science \
Večna pot 113, 1000 Ljubljana

#v(2.4em)

#align(center)[
    = Masters Thesis Topic Proposal
    *Noah Novšak*
]

#v(1.2em)

I, Noah Novšak, student of the 2nd cycle study program at the Faculty of Computer and Information Science, submit the following thesis topic proposal to be considered by the Committee for Student Affairs.

Working title:
#list(
    marker: "",
    indent: 1cm,
    tight: true,
    [Slovene: *Programski pristop k domnevi #PPT2*],
    [English: *A Software Approach to the #PPT2 Conjecture*]
)

This topic was already approved last year: *NO*

I declare that the mentors listed below have approved the submission of the thesis topic proposal described in the remainder of this document.

I would like to write the thesis in English with the following reason: I am a student of the Data Science program.


I propose the following mentor:
#list(
    tight: true,
    marker: "",
    indent: 1cm,
    [doc. dr. Aljaž Zalar],
    [Faculty of Computer and Information Science],
    [aljaz.zalar\@fri.uni-lj.si]
)

I propose the following co-mentor:
#list(
    tight: true,
    marker: "",
    indent: 1cm,
    [prof. dr. Igor Klep],
    [Faculty of Mathematics and Physics],
    [igor.klep\@fmf.uni-lj.si]
)

#align(right)[
    #align(bottom)[
        Ljubljana, December 5, 2024.
    ]
]

#pagebreak()

= Proposal of the thesis topic

== Narrow field of study

Quantum information theory, operator algebras, real algebraic geometry, mathematical optimization.

== Keywords

quantum communication,
positive-partial-transpose squared conjecture (#PPT2),
positive not completely positive maps,
semidefinite programming,
entanglement breaking,
separable states,
positive polynomials,
sums-of-squares polynomials

== Detailed thesis proposal

=== Introduction and problem formulation

We investigate the positive-partial-transpose squared (#PPT2) conjecture introduced by M. Christandl at Banff International Research Station Workshop: Operator Structures in Quantum Information Theory @Banff. The conjecture states that for any PPT map $Phi$, the composition $Phi compose Phi$ is entanglement breaking. Meaning, when the channel is applied to a part of an entangled state, it will result in a separable (non-entangled) state.

The importance of this conjecture comes to light in the context of quantum communication, for example when considering quantum repeaters. These seek to establish a secret key over a long distance out of entangled states over a smaller distance. While PPT states may allow the extraction of a private key @Horodecki_2005, the conjecture would imply that they cannot be used as a resource in a repeater @Christandl_2017.

In lower dimensions, the #PPT2 conjecture is known to hold. This is trivial for $n=2$ and proven twice independently for $n=3$ @Chen_2019 @Christandl_2019. In higher dimensions the conjecture has also been proven to hold for certain special cases, i.e. all Choi-type maps @Singh_2022 and Gaussian quantum channels @Christandl_2019. However, it is believed _not_ to hold in general, even though no proofs or counterexamples have been found.

The goal of this thesis is to investigate the #PPT2 conjecture in higher dimensions programmatically. Instead of attempting to find an analytic solution, we will generate random maps subject to the conjecture. Then we can reject the conjecture by counterexample or statistically validate it by testing separability. Testing separability is a well-studied problem in the context of polynomial sums-of-squares and semidefinite programming @Doherty_2004 @Fang_2020.

=== Related work

Besides the already mentioned work in lower dimensions and certain special cases, there have also been investigations in higher dimensions ($n gt.eq 4$). The authors of @jin2020 cover several possible approaches to finding a counterexample. The most direct approach is to find a PPT map $Phi$ and check if the composite map $Phi compose Phi$ is _not_ entanglement breaking. That is, finding the corresponding entanglement witness for the composite state. This is however extremely difficult. An alternative is to find a PPT entangled state as the composite channel and attempt to decompose it into two identical PPT channels. This saves us from verifying the composite is entanglement breaking, but in 4 dimensions still leaves us solving 256 nonlinear equations with 256 variables.

Determining separability is proven to be NP-hard. However, there are existing solutions @Doherty_2004 that provide separability criteria that can be cast as semidefinite programs (SDPs). This allows for a mathematically straightforward way to calculate a decomposition, given a strong enough computer. Although, it remains to be seen if this is feasible in higher dimensions.

The problem remains of finding the PPT map to test in the first place. Luckily, a method for generating positive maps that are not completely positive (PnCP) has already been proposed, based on the idea of generating positive polynomials that are not sums-of-squares @bhardwaj2020 @phdthesis. This is done by relaxing the positive condition to a multiplication with a polynomial, and then testing if the result is a sum of squares. This is a semidefinite program, and the main issue is choosing the right polynomial to make the SDP tractable. The proposed method is not yet scalable to higher dimensions due to the complexity of the underlying optimization problem and the performance issues of MATLAB.

=== Expected contributions

Our main goal is to provide a software library that can generate random PnCP maps in higher dimensions and use it to investigate the #PPT2 conjecture. Finding a counterexample would disprove the conjecture in higher dimensions, not finding one would still provide statistical evidence that if an example exists, it lies on the boundary of the cone of positive maps. Additionally, we will attempt to extend the library to generate maps that are $k$-positive and not ($k+1$)-positive in addition to PnCP maps.

=== Methodology

The existing `pncp` library will be rewritten in Julia, to include all the required tools that MATLAB provides, while increasing performance. To solve the necessary SDPs we will use the JuMP package, which provides an interface to performant solvers such as MOSEK and SeDuMi. The library will be tested against the existing MATLAB implementation to measure improvements in runtime, memory usage, and problem size scalability. If generating random maps in higher dimensions proves to be intractable, we will focus on optimizing the generation of specific subclasses of maps, such as ($k$, but not $k+1$)-positive maps.

#bibliography(
    style: "aps.csl",
    title: "References",
    "bibliography.bib"
)
