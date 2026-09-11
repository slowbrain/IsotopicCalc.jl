## IsotopicCalc

[![CI](https://github.com/slowbrain/IsotopicCalc.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/slowbrain/IsotopicCalc.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slowbrain/IsotopicCalc.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slowbrain/IsotopicCalc.jl)

Calculate isotopic patterns, exact monoisotopic masses, and candidate molecular formulas using bundled NIST isotope masses and abundances. Model ion transformations, isotope enrichment, and labeled mixtures; rank candidate formulas against measured isotope clusters and export the supporting evidence. Requires Julia 1.6 or later.

### Installation

```julia
using Pkg
Pkg.add("IsotopicCalc")
using IsotopicCalc
```

To use the code in a local checkout, run these commands from the repository directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=.
```

Then load the package with `using IsotopicCalc`. The registered release installed by `Pkg.add` may differ from this checkout.

### Quick start

```julia
using IsotopicCalc

mz = ion_mass("C3H6O"; ion="[M+H]+")
pattern = isotopic_pattern("C3H6O"; ion="[M+H]+", print=false)
candidates = find_formula(mz;
    ion="[M+H]+", tolerance_ppm=5,
    atom_pool=Dict("C"=>6, "H"=>14, "O"=>3), print=false)
```

| Task | API | Return value |
| --- | --- | --- |
| Neutral principal-isotope mass | `monoisotopic_mass` | `Float64` |
| Mass or m/z after atom and charge changes | `ion_mass` | `Float64` |
| Relative isotope pattern | `isotopic_pattern` | Vector of `(mass, abundance)` tuples |
| Pattern with isotope contributions and probability accounting | `isotopic_pattern(...; annotated=true)` | `PatternResult` |
| Candidate formulas | `find_formula` | Vector of `Compound` objects |
| Search with rejection counts and settings | `find_formula(...; diagnostics=true)` | `SearchResult` |
| Candidate ranking against a centroid cluster | `rank_formulas` | Vector of `RankedCandidate` objects |
| JSON-safe result snapshot | `result_dict` | Dictionary |

See [formula notation](#formula-notation), [ion transformations](#ion-transformations), [isotope labeling](#annotated-results-and-isotope-labeling), [search constraints](#search-constraints-and-diagnostics), and [cluster ranking](#ranking-measured-isotope-clusters) below.

### Compatibility notes

Version 0.8.0 uses explicit ion transformations and removes the legacy adduct API. Tuple output from `isotopic_pattern` and candidate-vector output from `find_formula` remain the defaults; annotation and diagnostics are opt-in.

Recent correctness fixes can change numerical results compared with older versions: monoisotopic masses use each element's principal isotope directly, isotope probabilities include permutation multiplicities, and finite resolving power produces merged peaks at abundance-weighted centroids. Use `R=Inf` to retain distinct isotope masses. Formula search also rejects impossible ion compositions and applies the corrected CHNO/halogen valence rules. See [CHANGELOG.md](CHANGELOG.md) for changes by version.

### Formula notation

These formulas all describe acetone:

```julia
"C3H6O"
"CH3COCH3"
"(CH3)2CO"
```

Parentheses may be nested; multipliers are applied directly to atom counts without expanding repeated text. Square brackets specify an isotope; `D` is shorthand for `[2H]`:

```julia
"[13C]C2H6O"  # One specified carbon-13 atom
"C3[2H]6O"   # Fully deuterated acetone
"C3D6O"      # Same composition
```

Formula and ion-notation arguments accept `AbstractString`, including `SubString` values. Formatted isotope formulas retain their square brackets. Malformed token boundaries and counts that overflow `Int` raise `ArgumentError`.

### Monoisotopic masses

```julia
monoisotopic_mass("C3H6O")             # ≈ 58.04186481295
monoisotopic_mass("Fe")                # 55.93493633 (56Fe)
monoisotopic_mass("[54Fe]")            # 53.93960899
monoisotopic_mass_protonated("C3H6O")  # ≈ 59.04914134518
```

Monoisotopic mass is the sum of the masses of the **most abundant isotope of each element**, honoring explicitly specified isotopes. It is computed directly from atom counts without generating a spectrum or rounding. It can differ from both the lightest peak and the most intense molecular peak of an isotope pattern.

### Isotopic patterns

```julia
isotopic_pattern("C3H6O")
```

At the default `R=10000`, the printed pattern is:

```text
Mass [amu]     Abundance [%]
58.04186481    100
59.04528960    3.351819
60.04653046    0.24406847
61.04959386    0.0069737949
```

By default, the return value is a `Vector{Tuple{Float64,Float64}}` containing unrounded mass (or m/z) and abundance relative to the strongest peak, whose abundance is `1.0`.

```julia
pattern = isotopic_pattern("C3H6O";
    abundance_cutoff=1e-5,
    R=10000,
    ion=IonSpec(),
    print=false,
    io=stdout,
)

# Resolve distinct isotope masses instead of merging nearby peaks:
fine_structure = isotopic_pattern("C3H6O"; R=Inf, print=false)

# Convenience function for a protonated pattern:
protonated = isotopic_pattern_protonated("C3H6O"; R=10000, print=false)
```

The calculation uses these rules:

- Isotope probabilities combine with their permutation multiplicities. For example, the probability of `k` carbon-13 atoms among `n` carbon atoms is `binomial(n,k) * p^k * (1-p)^(n-k)`.
- Resolving power is `R = m/Δm`. Starting with the lightest remaining peak, a group includes peaks within `first_mass/R` of that peak. Its mass is the abundance-weighted centroid, `sum(m*a)/sum(a)`, and its abundance is the sum of the constituent abundances. Groups do not grow by chaining adjacent spacings.
- The relative abundance cutoff is applied **after grouping and normalization**. Valid cutoffs are finite numbers in `[0,1]`; zero disables probability pruning, and one retains the strongest peak(s).
- For a positive cutoff, internal pruning discards at most `min(1e-12, abundance_cutoff*1e-8)` total absolute probability, apart from floating-point roundoff. Zero-cutoff calculations can become expensive for large formulas.
- Each element's tabulated abundance vector is normalized to sum to one. Explicit isotope labels have probability one.

### Ion transformations

Use the `ion` keyword in `ion_mass`, `isotopic_pattern`, and `find_formula`.
It accepts full `[M+...-...]z+/-` notation or an `IonSpec`. The signs **inside**
the brackets add or remove atoms; the suffix **outside** specifies charge.
Omit the suffix for an uncharged transformation. `M` denotes the input formula.

| Ion notation | Atom change | Charge |
|---|---|---:|
| `"[M]"` (default) | None | 0 |
| `"[M+H]+"` | Add one hydrogen (protonation) | +1 |
| `"[M-H]-"` | Remove one hydrogen (deprotonation) | −1 |
| `"[M-H]+"` | Remove one hydrogen, form a positive ion | +1 |
| `"[M-H]"` | Remove one hydrogen, remain uncharged | 0 |
| `"[M-2H]2-"` | Remove two hydrogens | −2 |
| `"[M+2H]2+"` | Add two hydrogens | +2 |
| `"[M+Na]+"` | Add sodium | +1 |
| `"[M+NH4]+"` | Add ammonium composition | +1 |
| `"[M+Cl]-"` | Add chlorine (chloride attachment) | −1 |
| `"[M-Cl]-"` | Remove chlorine | −1 |
| `"[M+H-H2O]+"` | Add hydrogen and lose water | +1 |
| `"[M]+"`, `"[M]2+"` | Remove electrons only | +1, +2 |
| `"[M]-"`, `"[M]2-"` | Add electrons only | −1, −2 |

```julia
ion_mass("C3H6O"; ion="[M-H]+")
isotopic_pattern("C3H6O"; ion="[M-H]-", print=false)
find_formula(ion_mass("C3H6O"; ion="[M-H]+"); ion="[M-H]+", print=false)

# The same hydrogen loss with an explicit, independent charge:
ion_mass("C3H6O"; ion=IonSpec(remove="H", charge=1))
ion_mass("C3H6O"; ion=IonSpec(remove="H", charge=-1))
ion_mass("C3H6O"; ion=IonSpec(remove="H"))
ion_mass("C3H6O"; ion=IonSpec(add="H", remove="H2O", charge=1))
ion_mass("C3H7O"; ion=IonSpec(charge=1)) # Formula already includes the ion's H
```

`IonSpec` stores only `additions`, `losses`, and signed `charge`. `string(spec)`
derives full notation from these fields; equivalent inputs may be normalized
(e.g. `[M+2H]2+` becomes `[M+H2]2+`). Atom changes are combined before validating
the resulting composition. Negative atom counts are rejected; isotope-labeled
atoms are separate entries and cannot silently replace unlabeled atoms.
Added unlabeled atoms use natural isotope distributions unless overridden.

For signed charge `z ≠ 0`, the conversion is
`m/z = (neutral_mass + atom_mass_change - z*electron_mass)/abs(z)`.
With `z = 0`, the result is `neutral_mass + atom_mass_change`.
`ion_mass` calculates this directly using reference principal-isotope masses;
isotope enrichment does not change those reference masses.

**Breaking API change in 0.8:** `adduct=` and compact strings such as `H+`, `H-`,
`Cl-`, `2H+2`, `+2`, and the empty string have been removed, with no aliases.
Use `ion="[M+H]+"`, `ion="[M-H]-"`, or the corresponding full notation.
For neutral calculations, omit `ion` or use `IonSpec()` / `"[M]"`.
`Compound.ion` replaces `Compound.adduct` and contains the full notation including
charge. Diagnostics use `:infeasible_ion`; JSON snapshots use schema version 2.

### Annotated results and isotope labeling

```julia
result = isotopic_pattern("C2";
    abundances=Dict("C"=>Dict("12C"=>0.8, "13C"=>0.2)),
    annotated=true, R=Inf, print=false)

[p.probability for p in result.peaks] # ≈ [0.64, 0.32, 0.04]
result.peaks[2].contributors[1].composition # Dict("12C"=>1, "13C"=>1)
result.retained_probability          # ≈ 1.0
```

With `annotated=true`, the return type is `PatternResult`. Each `IsotopePeak` has `mass`, relative `abundance`, absolute `probability`, and `contributors`. Each contributor records isotope counts, mass/m/z, absolute probability, and its population index. A merged peak's centroid and probability can be reconstructed from its contributors.

The result records retained probability, probability removed by the final cutoff, estimated probability removed by internal pruning, and the pruning bound. These account for the distribution up to floating-point roundoff. Settings include ion composition changes, source-population formulas/weights, abundance overrides, resolving power, cutoff, package version, and a SHA-256 fingerprint of the bundled isotope table. Annotation retains more information and uses more memory than tuple output.

Abundance overrides are **per call** and never modify the bundled table. Each element maps isotope labels to probabilities summing to one; omitted isotopes receive zero probability. Explicit isotope labels, including `D`, stay fixed. Enrichment changes pattern probabilities, while `monoisotopic_mass`, `ion_mass`, and formula search continue to use reference principal-isotope masses.

A mixture of molecular populations is different from independent per-atom enrichment:

```julia
populations = [
    IsotopePopulation("C2", 0.75; abundances=Dict("C"=>Dict("12C"=>1.0))),
    IsotopePopulation("[13C]2", 0.25),
]
result = isotopic_pattern("C2"; mixture=populations, annotated=true, print=false)
# Two peaks with absolute probabilities 0.75 and 0.25; no mixed 12C/13C molecule.
```

Mixture weights must sum to one. Population formulas must have equal elemental counts, ignoring isotope labels. Put abundance overrides inside each population when using a mixture. The base composition and every population receive the same ion specification.

### Finding formulas

```julia
matches = find_formula(59.0491;
    atom_pool=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10),
    tolerance_ppm=100,
    ion="[M+H]+",
    heuristics=true,
    print=false,
)
# Includes C3H6O [M+H]+ at m/z ≈ 59.049141, ppm ≈ -0.70
```

Atom-pool values are maximum counts. Results are `Compound` objects with `formula`, `ion`, `charge`, `mz`, and `ppm` fields, sorted by absolute ppm error. Formula strings use Hill ordering. Mass error is `1e6 * (observed_mz - calculated_mz)/calculated_mz`.

By default, for neutral closed-shell CHNO/halogen compositions, search requires integral, nonnegative `DBE = (2C + 2 + N - H - D - F - Cl - Br - I)/2`. Other element sets bypass this restricted valence model. Formula search is a candidate generator, not a general chemical-validity validator.

With `heuristics=true`, carbon-containing candidates also satisfy empirical limits `O ≤ 3C`, `N ≤ 4C`, `S ≤ C`, and `P ≤ 2C`. Set `heuristics=false` to search outside these limits, for example for carbon disulfide:

```julia
find_formula(75.9441423488;
    atom_pool=Dict("C"=>1, "S"=>2), tolerance_ppm=1, heuristics=false)
```

Search bounds each partial composition by the masses still reachable from its remaining atoms. Formula strings are constructed only for candidates within the requested mass-tolerance window. Larger atom pools and broader tolerance windows still increase search time. Inputs must have finite positive m/z, finite nonnegative ppm and absolute tolerances, and nonnegative atom counts.

`find_formula` and `isotopic_pattern` print by default. Set `print=false` for silent calls, or pass an `IO` destination such as `io=IOBuffer()` to capture formatted output. `isotopic_pattern_protonated` forwards pattern options, including `R`, `abundance_cutoff`, `annotated`, `abundances`, `mixture`, `print`, and `io`. Its ion specification is fixed to `ion="[M+H]+"`, so it rejects an additional `ion` argument. The direct mass functions do not print or accept spectrum settings.

### Search constraints and diagnostics

```julia
constraints = SearchConstraints(
    min_atoms=Dict("O"=>1),
    ratios=Dict(("H","C")=>(0, 3), ("S","C")=>(0, 2)),
    dbe=(0, 10),
)
report = find_formula(59.0491;
    ion="[M+H]+", tolerance_ppm=5, tolerance_da=0.001,
    constraints=constraints, diagnostics=true, print=false)
report.matches
report.diagnostics.status
report.diagnostics.rejections
```

The mass acceptance window uses the **larger** of the absolute tolerance and the ppm-derived tolerance: `abs(observed-calculated) ≤ max(tolerance_da, calculated*tolerance_ppm/1e6)`. `tolerance_da` is in mass/m/z units, so `0.001` is 1 mDa. Set either tolerance to zero to use only the other.

Explicit ratio intervals replace the corresponding default heuristic; a ratio with zero denominator is ignored. Minimum counts and explicit constraints apply even with `heuristics=false`. An explicit DBE constraint rejects compositions outside the supported CHNO/halogen valence model, rather than assigning them a guessed DBE.

`SearchResult` retains candidates and a snapshot of search settings. Diagnostic status is `:matched`, `:no_mass_match`, `:filtered_out`, or `:infeasible_ion`. `mass_matches` counts candidates inside the mass window **after atom-pool and atom-loss availability bounds**; rejected candidates are counted by their first failed constraint. It does not count every composition outside those bounds.

For a formula that already describes the ion, use `mode=:ion` with a charge-only specification:

```julia
find_formula(59.04914134518;
    mode=:ion, ion=IonSpec(charge=1),
    atom_pool=Dict("C"=>3, "H"=>7, "O"=>1), print=false)
# C3H7O, rather than the neutral precursor C3H6O
```

Ion mode's default CHNO/halogen model assumes even-electron ions, adding signed charge to the DBE numerator. For radical ions or compositions outside that assumption, set `SearchConstraints(valence=:none)`; this disables automatic parity/valence rejection, but any explicit DBE interval is still applied. Neutral mode searches precursor formulas and applies the ion transformation to each. Ion mode requires a charge-only specification because each candidate already describes the ion composition.

### Ranking measured isotope clusters

```julia
# Illustrative centroid mass/intensity pairs from one local cluster:
observed = [(44.02624, 1000.0), (45.0296, 22.1), (46.031, 2.2)]
precursor_mz = 44.02624
ranked = rank_formulas(precursor_mz, observed;
    search_options=(atom_pool=Dict("C"=>3,"H"=>8,"O"=>3), tolerance_ppm=1000),
    pattern_options=(R=1000,),
    match_tolerance_ppm=10, match_tolerance_da=0.05)

if !isempty(ranked)
    best = first(ranked)
    best.compound.formula
    best.isotope_similarity
    best.matches
end
```

`search_options` are forwarded to `find_formula`; `pattern_options` configure isotope prediction, including abundance overrides or mixtures. The search's ion specification is applied to each prediction. The ion, annotation, and printing options are controlled by the ranking function. Mixture formulas must be compatible with every candidate being ranked.

Use one local centroid cluster with finite nonnegative intensities, positive unique masses, and at least one nonzero intensity. Every candidate is evaluated against that same supplied cluster. Peak matching is ordered and one-to-one: maximize the number of pairs inside the matching window, then minimize summed squared normalized mass error. Missing predictions and unexplained measurements remain visible.

`RankedCandidate` returns:

- The `Compound`, precursor mass error, and a separate mass score.
- Cosine isotope-intensity similarity, including unmatched peaks in the vector norms.
- A least-squares intensity scale and per-pair mass/intensity residuals.
- Missing predicted-peak indices and unexplained observed-peak indices.
- The annotated prediction and a snapshot of measured data and matching settings.

Results sort by descending isotope similarity, then absolute precursor ppm error. The mass score is `exp(-0.5*(mass_error/window)^2)`, where the window is the search tolerance. Both scores describe agreement; neither is a calibrated probability of chemical identity. The matcher uses O(predicted peaks × observed peaks) memory, so supply a local cluster rather than an entire raw spectrum.

Observed-peak indices refer to the original input order; predicted-peak indices refer to `candidate.pattern.peaks`. Each `PeakMatch` exposes both indices, mass errors in Da and ppm, observed and fitted intensities, and the intensity residual. If formula search finds no candidates, ranking returns an empty vector.

### Saving results and running an example

```julia
using JSON
if !isempty(ranked)
    open("assignment.json", "w") do io
        JSON.print(io, result_dict(first(ranked)), 2)
    end
end
```

`result_dict` exports `PatternResult`, `SearchResult`, or `RankedCandidate` as a JSON-safe schema-version-2 snapshot. It includes settings and isotope-table fingerprint. Infinite resolving power or open constraint bounds are stored as `"Inf"`/`"-Inf"`. This is an export format; automatic deserialization into Julia result objects is not provided.

Run [examples/isotope_assignment.jl](examples/isotope_assignment.jl) for a complete synthetic-data example:

```sh
julia --project=. examples/isotope_assignment.jl assignment.json
```

### Development

Run the package tests from the repository directory:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

The suite covers isotope probabilities and resolution, formula parsing, ion atom/charge accounting, bounded formula search, enrichment and mixtures, diagnostic results, cluster matching, and JSON export. Package compatibility includes Julia 1.6 or later and JSON 0.21 or 1.x.

### Citing and references

See [CITATION.bib](CITATION.bib).

Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015), *Atomic Weights and Isotopic Compositions* (version 4.1). National Institute of Standards and Technology, Gaithersburg, MD. [NIST database](https://physics.nist.gov/Comp) (accessed August 22, 2023).
