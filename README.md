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

### Isotopic pattern distribution
Package uses isotopic abundances of elements from NIST online database (Atomic Weights and Isotopic Compositions with Relative Atomic Masses, see refenrence below).

#### Input formula
Acceptable sum formula formats are shown below (case for acetone):
```julia
C3H6O
CH3COCH3
(CH3)2CO
OH3(CH)3
```

There is also an option to specify concrete isotope using square bracket notation. 
```julia
[13C]H3COCH3    # acetone with one ¹³C isotope
C3[2H]6O        # fully deuterated acetone
C3D6O           # another notation for fully deuterated acetone
```
#### Examples
Main pupose of this package is to give isotopic pattern distribution for any chemical formula in following fashion:
```julia
julia> isotopicPattern("CH3CHOCH3");

Formula:  C3H6O
----------------------------
Mass [amu]     Abundance [%]
----------------------------
58.0419        100.0000
59.0452        3.2447
59.0461        0.0381
59.0481        0.0690
60.0461        0.2055
60.0486        0.0351
60.0494        0.0012
61.0495        0.0067
Found 8 isotopic masses for 1.0e-5 abundance limit.
```
Actual output of the function is a `Vector{Tuple{Float64,Float64}}` that contains molecular mass and its abundance relative to most abundant isotope. In this particular case for acetone output looks as follow:
```julia
8-element Vector{Tuple{Float64, Float64}}:
 (58.0419, 1.0)
 (59.0452, 0.0324471848781967)
 (59.0461, 0.0003809256493278667)
 (59.0481, 0.0006900793591262995)
 (60.0461, 0.002054993634531913)
 (60.0486, 0.0003509399355066257)
 (60.0494, 1.2359964968588418e-5)
 (61.0495, 6.667875838317436e-5)
 ```

The function has several keyword arguments used to customize the output. By default it has cut-off abundance set for 1e-5 as is always shown on function printout. The other keyword arguments include resolution (default value 10000) and possible adduct to the molecule (default set to "").
```julia
isotopicPattern(formula::String; 
                abundance_cutoff=1e-5, 
                R=4000, 
                adduct::String="", 
                print=true
               )
```
It's rather self explaining in following examples:
```julia
julia> isotopicPattern("CH3COCH3"; adduct = "Na+");

Formula:  C3H6O.Na+
----------------------------
Mass [amu]     Abundance [%]
----------------------------
81.0311        100.0000
82.0344        3.2447
82.0353        0.0381
82.0374        0.0690
83.0353        0.2055
83.0378        0.0351
83.0387        0.0012
84.0387        0.0067
Found 8 isotopic masses for 1.0e-5 abundance limit.

julia> isotopicPattern("CH3COCH3"; abundance_cutoff = 1e-8);

Formula:  C3H6O
----------------------------
Mass [amu]     Abundance [%]
----------------------------
58.0419        100.0000000
59.0452        3.2447185
59.0461        0.0380926
59.0481        0.0690079
60.0461        0.2054994
60.0486        0.0350940
60.0494        0.0012360
60.0515        0.0022391
60.0524        0.0000263
60.0544        0.0000198
61.0495        0.0066679
61.0519        0.0001265
61.0524        0.0001418
61.0528        0.0000134
61.0549        0.0000242
62.0528        0.0000721
62.0557        0.0000046
Found 17 isotopic masses for 1.0e-8 abundance limit.

julia> isotopicPattern("CH3COCH3"; R=4000);

Formula:  C3H6O
----------------------------
Mass [amu]     Abundance [%]
----------------------------
58.042         100.0000
59.045         3.2447
59.046         0.0381
59.048         0.0690
60.046         0.2055
60.049         0.0363
61.049         0.0067
Found 7 isotopic masses for 1.0e-5 abundance limit.
```
There are one special case implemented. Adduct `H⁺` is called inherently by `isotopicPatternProtonated`. 


### Monotisotopic mass calculation
For calculation of mono-istopic mass there is also `monoisotopicMass` and `monoisotopicMassProtonated`.
```julia
julia> monoisotopicMass("CH3COCH3")
58.0419
```

### Finding formulas
Rudimentary determination of compound sum formula from monoisotopic mass is implemented too.
```julia
julia> findFormula(58.0419);
Matching formulas:
 C3H6O  m/z: 58.041865  ppm: 0.61
 ```

The function has several keyword arguments used to customize the output. 
```julia
findFormula(mz_input::Float64;
            atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10),
            tolerance_ppm::Number=100,
            adduct::String="",
            charge::Int=0
           )
```
 By default algorithm searches within an interval of 100 ppm around monoisotopic mass and takes into consideration C, H, O, N to be possible building atoms. The associated numbers mean maximum amount of respective atoms to be considered.
So far adducts can be one of following `M+H`, `M+Na`, `M+K`, `M-H` and charge is to be expresed separately as in example.

```julia
julia> findFormula(46.000; tolerance_ppm=2000, charge=1);
Matching formulas:
 CH2O2+ m/z: 46.004931  ppm: -107.18
 NO2+   m/z: 45.992355  ppm: 166.23
 N2H2O+ m/z: 46.016164  ppm: -351.27
 CNH4O+ m/z: 46.02874   ppm: -624.4
 N3H4+  m/z: 46.039974  ppm: -868.24
 C2H6O+ m/z: 46.041316  ppm: -897.37
 CN2H6+ m/z: 46.05255   ppm: -1141.08
 C2NH8+ m/z: 46.065126  ppm: -1413.78
 C3H10+ m/z: 46.077702  ppm: -1686.32


julia> findFormula(59.0491; adduct="M+H", charge=1);
Matching formulas:
 C3H6O  [M+H]+  m/z: 59.049141  ppm: -0.7
 CN3H4  [M+H]+  m/z: 59.047799  ppm: 22.04
```


### Citing
See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

### References
Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015), Atomic Weights and Isotopic Compositions (version 4.1). [Online] Available: http://physics.nist.gov/Comp [2023, 08, 22]. National Institute of Standards and Technology, Gaithersburg, MD.

