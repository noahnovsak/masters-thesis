#let conf(
  title_en: "",
  title_sl: "",
  author: "",
  mentor: "",
  year: datetime.today().year(),
  abstract_en: [],
  abstract_sl: [],
  extended_abstract_sl: [],
  body,
) = {
  // Page layout according to UL FRI guidelines [cite: 109-118]
  set page(
    paper: "a4",
    margin: (
      top: 20mm,
      bottom: 30mm,
      inside: 30mm,
      outside: 20mm,
    ),
  )

  // Font settings: Titles in Helvetica/Arial, Text in Times New Roman [cite: 119]
  set text(font: "Libertinus Serif", size: 12pt, lang: "en")
  set par(leading: 0.8em, justify: true, linebreaks: "optimized") // Approx 1.3 line spacing [cite: 120]

  show heading: set text(font: "Helvetica", weight: "bold")
  show heading.where(level: 1): it => {
    pagebreak(weak: true, to: "odd") // Sections begin on odd pages [cite: 121]
    v(2em)
    it
    v(1em)
  }

  show math.equation: set text(font: "Libertinus Math")

  // --- 1. COVER PAGE (SLOVENE) [cite: 126-144] ---
  align(center)[
    #smallcaps("University of Ljubljana") \
    #smallcaps("Faculty of Computer and Information Science")

    #v(1fr)
    #text(size: 14pt, author) \
    #v(1em)
    #text(size: 16pt, weight: "bold", title_sl) \
    #v(1em)
    #text(size: 14pt, "MAGISTRSKO DELO") \
    #v(1em)
    #text(size: 12pt)[
      MAGISTRSKI ŠTUDIJSKI PROGRAM DRUGE STOPNJE \
      RAČUNALNIŠTVO IN INFORMATIKA \
      SMER: PODATKOVNE VEDE
    ]

    #v(1fr)
    #text(size: 12pt, "Mentor: " + mentor)

    #v(2em)
    #text(size: 12pt, "Ljubljana, " + str(year))
  ]
  pagebreak()
  pagebreak()

  // --- 2. TITLE PAGE (ENGLISH) [cite: 146, 157-164] ---
  align(center)[
    #smallcaps("University of Ljubljana") \
    #smallcaps("Faculty of Computer and Information Science")

    #v(1fr)
    #text(size: 14pt, author) \
    #v(1em)
    #text(size: 16pt, weight: "bold", title_en) \
    #v(1em)
    #text(size: 14pt, "MASTER'S THESIS") \
    #v(1em)
    #text(size: 12pt)[
      SECOND-CYCLE STUDY PROGRAMME \
      COMPUTER AND INFORMATION SCIENCE \
      TRACK: DATA SCIENCE
    ]

    #v(1fr)
    #text(size: 12pt, "Mentor: " + mentor)

    #v(2em)
    #text(size: 12pt, "Ljubljana, " + str(year))
  ]
  pagebreak()

  // --- 3. INTRODUCTORY PAGES (Unnumbered) [cite: 145] ---
  set page(numbering: none)

  // Intellectual Property Statement [cite: 148]
  heading(level: 1, numbering: none, outlined: false)[Intellectual Property Statement]
  [Statement text regarding intellectual property goes here...]
  pagebreak()

  // Acknowledgements [cite: 149]
  heading(level: 1, numbering: none, outlined: false)[Acknowledgements]
  [I would like to thank...]
  pagebreak()

  // English Abstract [cite: 152, 168]
  heading(level: 1, numbering: none, outlined: false)[Abstract]
  abstract_en
  pagebreak()

  // Slovenian Abstract [cite: 151, 169]
  heading(level: 1, numbering: none, outlined: false)[Povzetek]
  abstract_sl
  pagebreak()

  // Slovenian Extended Abstract (Min 10% of core length) [cite: 170-171]
  heading(level: 1, numbering: none, outlined: false)[Razširjeni povzetek]
  extended_abstract_sl
  pagebreak()

  // Table of Contents [cite: 153]
  outline(indent: auto)
  pagebreak()

  // --- 4. MAIN CONTENT (Numbered) [cite: 172] ---
  set page(numbering: "1", header: align(right)[#author: #title_en])
  counter(page).update(1)

  // Formatting for figures and tables [cite: 178-183]
  show figure: set align(center)
  show figure.caption: set text(size: 10pt)

  body
}
