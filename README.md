# A Software Approach to the PPT² Conjecture

Master's thesis by **Noah Novšak**, University of Ljubljana, Faculty of Computer
and Information Science (UL FRI).
Mentor: asist. prof. dr. Aljaž Zalar · Co-supervisor: prof. dr. Igor Klep.

The compiled thesis is [`thesis/main.pdf`](thesis/main.pdf).

## What this is

The **PPT² conjecture** (Christandl, 2012) asserts that the composition of any two
PPT maps is entanglement breaking. It is proven for maps on matrices up to
3 × 3 and for several structured families, but the general case is open; the
smallest open case, **4 × 4**, is the one this thesis attacks computationally.

This repository is a reproducible Julia pipeline that

1. **mass-produces entanglement witnesses** — provably indecomposable positive but
   not completely positive (PnCP) maps — via the Klep–McCullough–Šivic–Zalar
   construction, rationalizing every certificate after the SDP solve so each stored
   witness is *exact* (10,000 witnesses in under an hour);
2. **generates bound entangled PPT candidates** three ways — generic random
   sampling, partial-transpose-invariant sampling, and witness-guided extraction;
   and
3. **tests the conjecture** both by screening tens of thousands of composed channels
   with witness and DPS criteria, and by a see-saw SDP that searches the manifold of
   composed PPT maps directly.

No counterexample is found. The central finding: every one of the 10,000 witnesses
attains a negative optimum somewhere on the PPT cone, yet not one fires on the
composition manifold — suggesting the conjecture holds in dimension this case.

## Repository layout

```
code/            Julia package `ppt2` + the search pipeline
  src/           Library: PnCP construction (pncp.jl), PPT/PPT² SDPs (ppt2.jl), I/O (io.jl)
  scripts/       Command-line drivers for the long-running jobs — see code/scripts/README.md
  notebooks/     Interactive workflows (rationalization, PPT-state and UPB sampling, …)
  test/          Test suite: construction reference values, positivity checks, scripts
  results/       The seed-42 final scan: outputs (.jld2 / .csv), logs, and run_all.sh
thesis/          Typst source (main.typ), bibliography, figures, and compiled main.pdf
  template/      Two interchangeable layouts (paper / book) + citation styles
LICENSE.md       CC BY-SA 4.0 for the text; GPL-3.0 for the code and results
```

## Running the code

Everything runs from the `code/` directory (the Julia project root). The SDP steps
use **Mosek**, so a valid Mosek license is required, but other solvers (e.g., **Hypatia**) can be substituted.

```sh
cd code
# install dependencies once
julia --project=. -e 'using Pkg; Pkg.instantiate()'
# run the test suite
julia --project=. -e 'using Pkg; Pkg.test()'
```

The generation and testing scripts are multithreaded, set the thread count with
`-t`:

```sh
julia --project=. -t auto scripts/gen_pncp.jl --total 10000 --batch 200 -n 4 -m 4
```

See [`code/scripts/README.md`](code/scripts/README.md) for the full pipeline:
each driver, its options, the resume/reproducibility model, and the output formats.

### Reproducing the thesis results

[`code/results/run_all.sh`](code/results/run_all.sh) is the end-to-end driver that
produced the figures and tables in the thesis (seed 42, 4 × 4). It runs the steps
sequentially, logging each to `results/logs/` and its wall time to a timings table,
and is resumable. Completed steps are skipped on a rerun.

```sh
cd code
bash results/run_all.sh
```

The committed `code/results/` directory already holds that scan's outputs; the
`logs/` and `driver.out` files are the original run records.

## Building the thesis

The thesis is written in [Typst](https://typst.app). From the repository root:

```sh
typst compile thesis/main.typ
```

This produces `thesis/main.pdf`. The two layouts in `thesis/template/` are
interchangeable, switch the import at the top of `main.typ` between `book.typ` and
`paper.typ`. (If Typst warns about a missing math font, install the document fonts
it names, or let it fall back.)

## License

The text, figures, and results are released under
**Creative Commons Attribution-ShareAlike 4.0**; the source code and the software
developed for the thesis are released under the **GNU General Public License,
version 3 (or newer)**. See [LICENSE.md](LICENSE.md) for details. Please cite the
thesis and its author when reusing this work.
