# IsotopicCalc

<!-- [![Coverage](https://codecov.io/gh/slowbrain/IsotopicCalc.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slowbrain/IsotopicCalc.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/I/IsotopicCalc.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html) -->

This package is first relatively simple atempt to calculate monoisotopic mass and/or isotopic pattern for any given chemical formula. The reverse function to determine the formula from given molecular mass is yet to be implemented.

To install run one of the following commands:
```julia
]add IsotopicCalc
```
or
```julia
using Pkg; Pkg.add("IsotopicCalc.jl")
```

## Usage
```julia
using IsotopicCalc

julia> monoisotopicMass("C3H6O")
58.041865

julia> isotopicPattern("CH3COCH3");

Formula:  C3H6O
-------------------------------
Mass [amu]      Abbundance [%]
-------------------------------
58.041865       100.0
59.04522        3.244718
60.04611        0.205499
59.048142       0.069008
59.046082       0.038093
60.048575       0.035094
61.049465       0.006668
60.049437       0.001236
Found 8 isotopic masses for 1.0e-5 abundance limit.
```

### Citing
See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
