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

This topic was already approved last year: *YES*

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
        Ljubljana, December 5, 2025.
    ]
]

#pagebreak()

= Key-words
quantum communication,
positive-partial-transpose squared conjecture (#PPT2),
positive not completely positive maps,
semidefinite programming,
entanglement breaking,
separable states,
positive polynomials,
sums-of-squares polynomials

= Detailed thesis proposal
*Problem & State of the Art* The #PPT2 conjecture posits that the composition of any positive-partial-transpose (PPT) map with itself is entanglement breaking. This hypothesis is critical for quantum communication, particularly regarding the security and feasibility of quantum repeaters @Horodecki_2005 @Christandl_2017. While the conjecture is proven for low dimensions ($n <= 3$) @Chen_2019 @Christandl_2019 and specific classes such as Choi type maps @Singh_2022 and Gaussian quantum channels @Christandl_2019, it is widely believed to fail in higher dimensions ($n >= 4$). However, no general proof or counterexample exists. A MATLAB-based library exists for generating map candidates based on positive polynomials that are not sums-of-squares @Bhardwaj_2023 @Bhardwaj_2020, but it is limited in scalability and performance due to the complexity of the underlying semidefinite programs (SDPs) and MATLAB's inefficiencies.

*Expected Contributions / Technical Outcome* This work will deliver a high-performance, open-source software library in Julia capable of generating random positive maps that are not completely positive (PnCP) in higher dimensions. The primary contribution is a rigorous computational stress-test of the #PPT2 conjecture: we aim to either isolate a numerical counterexample, thereby disproving the conjecture, or provide significant statistical evidence of its validity. Additionally, the project will expand the tooling landscape by implementing methods to generate maps that are strictly $k$-positive but not $(k+1)$-positive, providing new resources for analyzing the geometry of quantum channels.

*Methodology & Validation* We will re-architect the existing generation framework using Julia and the JuMP package to interface with high-performance solvers like MOSEK. The approach relies on polynomial sum-of-squares (SOS) relaxations to generate map candidates, followed by semidefinite programming (SDP) to verify separability criteria @Doherty_2004 @Harrow_2017: first, we will benchmark the new library against the legacy MATLAB implementation to quantify improvements in runtime, memory usage, and dimensional scalability. Second, we will deploy large-scale randomized sampling to probe the boundaries of the positive map cone, using SDP certificates to validate the entanglement-breaking properties of the composite maps.


#bibliography(
    style: "aps.csl",
    title: "References",
    "bibliography.bib"
)
