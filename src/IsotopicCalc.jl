"""
    IsotopicCalc

A Julia package for calculating isotopic patterns and monoisotopic masses of chemical compounds.

# Features
- Calculate complete isotopic pattern distributions with natural abundances
- Compute monoisotopic masses for any chemical formula
- Find molecular formulas from experimental m/z values
- Support for complex formula notations (parentheses, specific isotopes)
- Handle various ionization states and adducts (H⁺, Na⁺, K⁺, etc.)
- Uses NIST atomic weights and isotopic composition data

# Main Functions
- [`isotopicPattern`](@ref): Calculate full isotopic distribution
- [`monoisotopicMass`](@ref): Calculate monoisotopic mass
- [`findFormula`](@ref): Find formulas matching an m/z value
- [`isotopicPatternProtonated`](@ref): Convenience function for [M+H]⁺
- [`monoisotopicMassProtonated`](@ref): Monoisotopic mass of [M+H]⁺
- [`Compound`](@ref): Data structure for formula search results

# Example Usage
```julia
using IsotopicCalc

# Calculate isotopic pattern for acetone
isotopicPattern("CH3COCH3")

# Get monoisotopic mass
mass = monoisotopicMass("C6H12O6")  # Returns 180.0634

# Find formulas matching an m/z value
matches = findFormula(58.0419; tolerance_ppm=10)

# Calculate pattern with adduct
isotopicPattern("C6H12O6"; adduct="Na+")
```

# Formula Notation
Supported formats:
- Simple: `"C3H6O"`, `"H2SO4"`
- Parentheses: `"(CH3)2CO"`, `"Ca(OH)2"`
- Specific isotopes: `"[13C]H3COCH3"`, `"C3[2H]6O"`
- Deuterium shorthand: `"C3D6O"` (same as `"C3[2H]6O"`)

# References
Based on NIST Atomic Weights and Isotopic Compositions database (version 4.1).
Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015),
http://physics.nist.gov/Comp

See also: README.md for detailed examples and CITATION.bib for citing this package.
"""
module IsotopicCalc

    using JSON, Printf

    const ELEMENTS = JSON.parsefile(Base.Filesystem.joinpath(dirname(@__FILE__),"elements.json"));

    include("isotopicPattern.jl")
    include("findFormula.jl")

    export isotopicPattern, monoisotopicMass, isotopicPatternProtonated, monoisotopicMassProtonated, findFormula, Compound

end # module
