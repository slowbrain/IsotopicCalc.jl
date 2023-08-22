## IsotopicCalc

<!-- [![Coverage](https://codecov.io/gh/slowbrain/IsotopicCalc.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slowbrain/IsotopicCalc.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/I/IsotopicCalc.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html) -->

### Installation
To install run one of the following commands:
```julia
]add IsotopicCalc
```
or
```julia
using Pkg; Pkg.add("IsotopicCalc.jl")
```

### Usage
Below are examples of 

```julia
using IsotopicCalc

julia> monoisotopicMass("C3H6O")
58.041865

julia> isotopicPattern("CH3COCH3");

Formula:  CH3COCH3.H+
-------------------------------
Mass [amu]      Abundance [%]
-------------------------------
59.0491         100.000
60.0525           3.245
60.0534           0.038
60.0554           0.081
61.0534           0.205
61.0559           0.035
61.0567           0.001
62.0567           0.007
Found 8 isotopic masses for 1.0e-5 abundance limit.
```

### Citing
See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

