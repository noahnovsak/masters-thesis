#let conf(
  title_en: "",
  title_sl: "",
  author: "",
  mentor: "",
  cosupervisor: none,
  year: datetime.today().year(),
  keywords_en: "",
  keywords_sl: "",
  abstract_en: [],
  abstract_sl: [],
  extended_abstract_sl: [],
  acronyms: none,
  acknowledgements: [I would like to thank...],
  dedication: none,
  code_url: none,
  body,
) = {
  // Page layout per UL FRI guidelines: A4, bound edition (inner margin wider).
  set page(
    paper: "a4",
    margin: (top: 20mm, bottom: 30mm, inside: 30mm, outside: 20mm),
    numbering: none,
  )

  // Headings in Helvetica/Arial, body in a Times-like serif, ~1.3 line spacing.
  set text(font: "Times New Roman", size: 12pt, lang: "en")
  set par(leading: 0.8em, justify: true, linebreaks: "optimized")

  show heading: set text(font: "Helvetica", weight: "bold")
  show heading.where(level: 1): it => {
    pagebreak(weak: true, to: "odd") // Chapters begin on odd pages.
    v(2em)
    it
    v(1em)
  }
  set heading(numbering: "1.1")

  set math.equation(numbering: "(1)")
  show math.equation: set text(font: "New Computer Modern Math")

  // Centered figures with bold caption labels ("Figure N", "Table N").
  show figure: set align(center)
  show figure.caption: it => context {
    set text(size: 10pt)
    strong[#it.supplement #it.counter.display(it.numbering)#it.separator]
    it.body
  }

  // --- COVER PAGE (Slovene) ---
  {
    set align(center)
    text(size: 14pt, smallcaps[Univerza v Ljubljani]); linebreak()
    text(size: 14pt, smallcaps[Fakulteta za računalništvo in informatiko])
    v(9em)
    text(size: 14pt, author); linebreak()
    v(1em)
    text(size: 18pt, weight: "bold", title_sl)
    v(1.2em)
    text(size: 14pt, weight: "bold")[MAGISTRSKO DELO]
    v(0.4em)
    text(size: 12pt)[
      MAGISTRSKI ŠTUDIJSKI PROGRAM DRUGE STOPNJE \
      RAČUNALNIŠTVO IN INFORMATIKA \
      SMER: PODATKOVNE VEDE
    ]
    v(1fr)
    text(size: 12pt)[Mentor: #mentor]
    if cosupervisor != none {
      linebreak()
      text(size: 12pt)[Somentor: #cosupervisor]
    }
    v(1.5em)
    text(size: 12pt)[Ljubljana, #str(year)]
  }
  pagebreak()

  // --- TITLE PAGE (English) ---
  {
    set align(center)
    text(size: 14pt, smallcaps[University of Ljubljana]); linebreak()
    text(size: 14pt, smallcaps[Faculty of Computer and Information Science])
    v(9em)
    text(size: 14pt, author); linebreak()
    v(1em)
    text(size: 18pt, weight: "bold", title_en)
    v(1.2em)
    text(size: 14pt, weight: "bold")[MASTER'S THESIS]
    v(0.4em)
    text(size: 12pt)[
      THE 2ND CYCLE MASTER'S STUDY PROGRAMME \
      COMPUTER AND INFORMATION SCIENCE \
      TRACK: DATA SCIENCE
    ]
    v(1fr)
    text(size: 12pt)[Supervisor: #mentor]
    if cosupervisor != none {
      linebreak()
      text(size: 12pt)[Co-supervisor: #cosupervisor]
    }
    v(1.5em)
    text(size: 12pt)[Ljubljana, #str(year)]
  }
  pagebreak()

  // --- COPYRIGHT / LICENSE PAGE (CC BY-SA 4.0) ---
  {
    set text(size: 10pt)
    v(1fr)
    par(justify: true)[
      This work is licensed under the _Creative Commons Attribution–ShareAlike
      4.0 International (CC BY-SA 4.0)_ license. This means that the text, images,
      graphs, and other components of this work may be freely shared, reproduced,
      made available to the public, and adapted for any purpose, including
      commercial use, provided that the author is clearly credited (preferably also
      the title and a link to the original), a link to the license is provided, any
      changes are indicated, and adaptations are distributed under the same license
      (CC BY-SA 4.0). The license does not permit the addition of extra legal or
      technological restrictions and does not apply to parts for which the rights
      holder is not the author. Details of the license are available at
      #link("https://creativecommons.org").
    ]
    v(0.5em)
    align(center, grid(
      columns: 3,
      gutter: 1em,
      image("license/cc_cc_30.svg", height: 0.9cm),
      image("license/cc_by_30.svg", height: 0.9cm),
      image("license/cc_sa_30.svg", height: 0.9cm),
    ))
    if code_url != none {
      v(1.2em)
      par(justify: true)[
        The source code of the thesis, its results, and the software developed for
        this purpose are licensed under the GNU General Public License, version 3
        (or later). This means it may be freely distributed and/or modified under
        its terms, and is available at #link(code_url). Details of the license are
        available at #link("https://www.gnu.org/licenses/").
      ]
    }
    v(1.5em)
    align(center, text(size: 9pt, smallcaps[© #str(year) #author]))
  }
  pagebreak()

  // --- ACKNOWLEDGEMENTS ---
  {
    align(center, text(size: 16pt, weight: "bold", smallcaps[Acknowledgements]))
    v(0.8cm)
    emph(acknowledgements)
    v(0.8cm)
    align(right, emph[#author, #str(year)])
  }
  pagebreak()

  // --- DEDICATION (optional) ---
  if dedication != none {
    v(20%)
    align(right, block(width: 60%, emph(dedication)))
    pagebreak()
  }

  // --- TABLE OF CONTENTS ---
  outline(title: [Contents], depth: 2, indent: auto)
  pagebreak()

  // --- LIST OF ACRONYMS (optional) ---
  if acronyms != none {
    heading(level: 1, numbering: none, outlined: false)[List of acronyms]
    acronyms
  }

  // --- ABSTRACT (English) ---
  heading(level: 1, numbering: none, outlined: false)[Abstract]
  [*Title:* #title_en]
  v(0.5em)
  abstract_en
  v(1em)
  text(weight: "bold")[Keywords]
  linebreak()
  emph(keywords_en)

  // --- POVZETEK (Slovene) ---
  heading(level: 1, numbering: none, outlined: false)[Povzetek]
  [*Naslov:* #title_sl]
  v(0.5em)
  abstract_sl
  v(1em)
  text(weight: "bold")[Ključne besede]
  linebreak()
  emph(keywords_sl)

  // --- EXTENDED ABSTRACT (Slovene, required for theses written in English) ---
  heading(level: 1, numbering: none, outlined: false)[Razširjeni povzetek]
  extended_abstract_sl

  // --- MAIN CONTENT (numbered, with running headers) ---
  counter(page).update(1)
  set page(
    numbering: "1",
    footer: none,
    header: context {
      let p = here().page()
      let chaps = query(heading.where(level: 1))
        .filter(h => h.numbering != none and h.location().page() <= p)
      let chap = if chaps.len() > 0 { upper(chaps.last().body) } else { [] }
      let pnum = counter(page).display("1")
      set text(size: 9pt, style: "italic")
      block(below: 4pt, width: 100%, if calc.even(p) {
        grid(columns: (auto, 1fr), pnum, align(right, chap))
      } else {
        grid(columns: (1fr, auto), chap, align(right, pnum))
      })
      line(length: 100%, stroke: 0.5pt)
    },
  )

  body
}
