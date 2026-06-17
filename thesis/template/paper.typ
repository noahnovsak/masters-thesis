// paper.typ — Typst port of the UL FRI "friteza" thesis-paper template
// (style/friteza.cls, based on PNAS pnas-new.cls). Reproduces the
// conference/journal-paper layout: single-column body, single-column title
// block and abstract, a left-margin sidebar carrying the supervisor/copyright
// metadata and the other-language title/abstract/keywords, an alternating
// running footer, and PNAS-style headings.
//
// Mirrors the role of user_defs.tex + the class file; write the thesis in
// main.typ. The book-style `book.typ` is kept for reference.

#import "@preview/drafting:0.2.2": *
#import "@preview/algorithmic:1.0.7": *
#import "@preview/lemmify:0.1.8": *

// --- Colors lifted from friteza.cls / pnasresearcharticle.sty ---
#let col-accent  = rgb(240, 14, 14)    // P1797 — author + keywords
#let col-title   = rgb(51, 51, 51)     // black80 — title text
#let col-rule    = rgb(128, 128, 128)  // black50 — hrules

// Fonts matching the LaTeX template exactly:
//  - serif: lmodern  -> Latin Modern Roman
//  - sans:  helvet   -> Helvetica (scaled to 0.95, helvet's [scaled] default)
//  - math:  lmodern  -> Latin Modern Math
#let serif = "Latin Modern Roman"
#let sans = "Helvetica"
#let mathfont = "Latin Modern Math"

// helvet [scaled] downscales Helvetica so its x-height matches the serif.
#let hscale = 0.95
#let sf(s) = s * hscale  // scaled sans size

// Theorem/definition style: the label (type, number, name) is set in sans
// (Helvetica), echoing the run-in headings; the body stays in the serif body
// font. Mirrors lemmify's `thm-style-simple` but swaps the label font.
#let thm-style-heading(thm-type, name, number, body) = block(width: 100%, breakable: true)[#{
  set align(left)
  text(font: sans, weight: "bold")[#thm-type#if number != none [ #number]]
  if name != none {
    text(font: serif, style: "italic")[ (#name)]
  }
  h(0.5em)
  body
}]

#let (theorem, definition, rules: thm-rules) = default-theorems(
  "thm-group", lang: "en",
  thm-numbering: thm-numbering-heading.with(max-heading-level: 1),
  max-reset-level: 1,
  thm-styling: thm-style-heading,
)

// Keywords are stored comma-separated; the template shows them pipe-separated.
#let fmt-keywords(kw) = kw.split(", ").join(" | ")

// Bibliography renderer (called from main.typ): two-column PNAS layout.
// Path is relative to this file, so the .bib sits one directory up.
#let render-bib(path: "../bibliography.bib") = {
  linebreak()
  columns(2, {
    heading(level: 1, numbering: none)[Bibliography]
    bibliography(title: none, style: "pnas.csl", path)
  })
}

#let conf(
  title_en: "",
  title_sl: "",
  author: "",
  mentor: "",
  cosupervisor: none,
  program_en: "Master's study programme Computer and Information Science",
  program_sl: "Magistrski študijski program Računalništvo in informatika",
  track_en: "Data Science",
  track_sl: "Podatkovne vede",
  year: datetime.today().year(),
  keywords_en: "",
  keywords_sl: "",
  abstract_en: [],
  abstract_sl: [],
  extended_abstract_sl: [],
  code_url: none,
  acknowledgements: none,
  body,
) = {
  // Sidebar/gutter geometry: how far left-margin content shifts out of the body
  // and how wide it is. The footer stretches by the same shift so its right edge
  // still lands on the body's right margin. The logo hangs slightly further out.
  let gutter-shift = 157pt
  let gutter-width = 145pt

  // --- Page geometry: wide left margin holds the sidebar (PNAS 2.25in). ---
  set page(
    paper: "a4",
    margin: (left: 70mm, right: 15mm, top: 16mm, bottom: 15mm),
    numbering: "1",
    // Running footer; alternates between recto (odd) and verso (even) pages
    // and stretches across the full page width (into the left margin), as the
    // LaTeX template does with \hspace*{-2.25in}.
    footer: context {
      set text(font: sans, size: sf(7pt))
      let sep = h(6pt) + sym.bar.v + h(6pt)
      // Author rendered "Lastname F." for the running footer (last token is the
      // surname; every preceding given name is reduced to an initial).
      let parts = author.split(" ")
      let lead = emph[#parts.last() #parts.slice(0, -1).map(p => p.clusters().first() + ".").join(" ")]
      if here().page() == 1 {
        // First page: the full page range stands in for the page number, and
        // the author lead is omitted.
        let range = text(weight: "bold", [1#sym.dash.en#counter(page).final().first()])
        let meta = [BMA-RI-PV#sep Master's thesis#sep #year#sep #range]
        move(dx: -gutter-shift, box(width: 100% + gutter-shift, align(right, meta)))
      } else {
        let pagenum = text(weight: "bold", counter(page).display("1"))
        let meta = [BMA-RI-PV#sep Master's thesis#sep #year#sep #pagenum]
        // odd (recto): author | … meta+page (outer-right);  even (verso): page | … author
        let (lft, rgt) = if calc.odd(here().page()) { (lead, meta) } else { (pagenum, lead) }
        move(dx: -gutter-shift, box(width: 100% + gutter-shift, grid(
          columns: (1fr, auto),
          align(left, lft),
          align(right, rgt),
        )))
      }
    },
  )

  // --- Body text: 9pt serif, justified, dense leading (documentclass[9pt]). ---
  set text(font: serif, size: 9pt, lang: "en")
  set par(leading: 0.62em, justify: true, first-line-indent: 1em, spacing: 0.62em)

  // --- Headings: sans-serif, PNAS numbering, run-in for deeper levels. ---
  set heading(numbering: "1.1")
  show heading: set text(font: sans)

  show heading.where(level: 1): it => {
    set text(size: sf(11pt), weight: "bold")
    block(above: 1.6em, below: 0.9em)[
      #if it.numbering != none [#counter(heading).display("1").#h(0.4em)]
      #it.body
    ]
  }
  // Levels 2-4 are run-in: the heading sits inline and the paragraph
  // continues on the same line (returning inline content, no block).
  show heading.where(level: 2): it => {
    v(1.4em, weak: true)
    set text(size: sf(9.5pt), weight: "bold")
    if it.numbering != none [#numbering("1.1", ..counter(heading).at(it.location())).#h(0.4em)]
    it.body
    [.]
    h(0.6em)
  }
  show heading.where(level: 3): it => {
    v(1.4em, weak: true)
    set text(size: sf(9pt), weight: "bold", style: "italic")
    if it.numbering != none [#numbering("1.1.1", ..counter(heading).at(it.location())).#h(0.4em)]
    it.body
    [.]
    h(0.6em)
  }
  show heading.where(level: 4): it => {
    set text(size: sf(9pt), weight: "bold", style: "italic")
    it.body
    [.]
    h(0.6em)
  }

  // --- Math, code, links ---
  set math.equation(numbering: "(1)")
  show math.equation: set text(font: mathfont)
  // Latin Modern Mono is not bundled with Typst; New Computer Modern Mono is its successor.
  show raw: set text(size: 8.5pt, font: "Latin Modern Mono")

  // --- References: hierarchical section number for headings (the prose
  // supplies the word "section"), equation number for equations. ---
  show ref: it => {
    let el = it.element
    if el == none {
      it
    } else if el.func() == math.equation {
      let num = numbering(el.numbering, ..counter(math.equation).at(el.location()))
      link(el.location())[#num]
    } else {
      it
    }
  }

  // --- References: sans-serif, 6.5pt, tight, matching the PNAS class. ---
  show bibliography: set text(font: sans, size: sf(6.5pt))
  show bibliography: set par(leading: 0.5em, spacing: 0.5em)

  // --- Figures / tables: sans, small, bold "Figure N." label. ---
  show figure: set align(center)
  // Extra breathing room below a figure before the body text resumes. Scoped to
  // real figures (image / table / algorithm); lemmify renders theorems and
  // definitions as figures too (kind "thm-group"), which must keep tight spacing.
  show figure.where(kind: image): set block(below: 1.8em)
  show figure.where(kind: table): set block(below: 1.8em)
  show figure.where(kind: "algorithm"): set block(below: 1.8em)
  show figure.caption: it => context {
    set text(font: sans, size: sf(8pt))
    set par(first-line-indent: 0pt)
    [#strong[#it.supplement #it.counter.display(it.numbering)#it.separator]#it.body]
  }

  // --- Sidebar metadata, reused below: faculty / supervisor / copyright. ---
  let metadata = block(width: gutter-width)[
    #set text(font: sans, size: sf(7pt), fill: col-title)
    #set par(first-line-indent: 0pt, justify: true, leading: 0.5em)
    *Faculty of Computer and Information Science* \
    #v(4pt)
    *Supervisor:* #mentor \
    #if cosupervisor != none [*Co-supervisor:* #cosupervisor \ ]
    #v(8pt)
    #text(size: sf(7pt))[
      *Copyright:* This work is licensed under a CC BY-SA 4.0 license.
      #if code_url != none [ The software and results produced in this work are released under GPL-3.0-or-later at #link(code_url).]
    ]
  ]

  // ============================ FIRST-PAGE SIDEBAR ============================
  // Logo at the very top of the wide left margin; hangs slightly further out
  // than the text sidebar.
  place(top + left, dx: -194pt, image("figures/UL-FRI.png", width: 84pt))

  // Slovene title/abstract/keywords, plain, aligned to the page bottom.
  place(bottom + left, dx: -gutter-shift, dy: 0mm, block(width: gutter-width)[
    #set text(font: sans, size: sf(7pt), fill: col-title)
    #set par(first-line-indent: 0pt, justify: true, leading: 0.5em)
    #strong(title_sl)
    #v(4pt)
    #abstract_sl
    #v(5pt)
    #fmt-keywords(keywords_sl)
  ])

  // ============================ TITLE BLOCK ============================
  {
    set par(first-line-indent: 0pt, justify: false, leading: 0.7em)
    text(font: sans, size: sf(20pt), weight: "bold", fill: col-title, title_en)
    v(8pt)
    text(font: sans, size: sf(9pt), weight: "bold", fill: col-accent, author)
    v(4pt)
    text(font: sans, size: sf(7pt), fill: col-title)[#program_en, study field #track_en]
    v(10pt)
  }

  // Place the metadata in the margin level with the abstract: the current flow
  // position is now at the abstract's top; subtract the top page margin to get
  // the container-relative offset for a top-anchored placement.
  context place(top + left, dx: -gutter-shift, dy: here().position().y - 16mm, metadata)

  // ============================ ABSTRACT (full width, sans bold) ============================
  {
    set text(font: sans, size: sf(8pt), weight: "bold", fill: col-title)
    set par(first-line-indent: 0pt, justify: true, leading: 0.6em)
    abstract_en
    v(6pt)
    text(size: sf(7pt), weight: "regular", fill: col-accent, fmt-keywords(keywords_en))
  }
  v(12pt)

  // ============================ BODY (single column) ============================
  // The main body stays single column; only the references are set in two
  // columns (handled at the bibliography call in main.typ).
  // The extended Slovene abstract is intentionally omitted in this format;
  // `extended_abstract_sl` is accepted for compatibility but not rendered.
  body
}
