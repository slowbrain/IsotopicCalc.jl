"""
    IsotopicCalc

A Julia package for calculating isotopic patterns and monoisotopic masses of chemical compounds.

# Features
- Calculate complete isotopic pattern distributions with natural abundances
- Compute monoisotopic masses for any chemical formula
- Find molecular formulas from experimental m/z values
- Support for complex formula notations (parentheses, specific isotopes)
- Model atom additions, atom losses, and signed charge with explicit ion transformations
- Uses NIST atomic weights and isotopic composition data

# Main Functions
- [`isotopic_pattern`](@ref): Calculate full isotopic distribution
- [`monoisotopic_mass`](@ref): Calculate monoisotopic mass
- [`find_formula`](@ref): Find formulas matching an m/z value
- [`isotopic_pattern_protonated`](@ref): Convenience function for [M+H]⁺
- [`monoisotopic_mass_protonated`](@ref): Monoisotopic mass of [M+H]⁺
- [`Compound`](@ref): Data structure for formula search results
- [`IonSpec`](@ref), [`ion_mass`](@ref): Explicit ion compositions and direct m/z
- [`PatternResult`](@ref), [`IsotopePopulation`](@ref): Annotated patterns and labeling mixtures
- [`SearchConstraints`](@ref), [`SearchResult`](@ref): Constrained searches and diagnostics
- [`rank_formulas`](@ref): Candidate ranking against observed isotope clusters
- [`result_dict`](@ref): JSON-safe reproducibility snapshots

# Example Usage
```julia
using IsotopicCalc

# Calculate isotopic pattern for acetone
isotopic_pattern("CH3COCH3")

# Get monoisotopic mass
mass = monoisotopic_mass("C6H12O6")  # Returns 180.0634

# Find formulas matching an m/z value
matches = find_formula(58.0419; tolerance_ppm=10)

# Calculate the pattern for sodium attachment
isotopic_pattern("C6H12O6"; ion="[M+Na]+")
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

    using JSON, SHA, TOML

    const PACKAGE_VERSION = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
    const ISOTOPE_DATA_SHA256 = bytes2hex(sha256(read(joinpath(@__DIR__, "elements.json"))))

    # Load element data that will be used by submodules
    const ELEMENTS = JSON.parsefile(Base.Filesystem.joinpath(dirname(@__FILE__),"elements.json"));

    # Include and load submodules
    include("IsotopicPattern.jl")
    include("FindFormula.jl")

    # Import and re-export functions from submodules
    using .IsotopicPatternModule: isotopic_pattern, monoisotopic_mass, isotopic_pattern_protonated, monoisotopic_mass_protonated
    using .FindFormulaModule: Compound, find_formula, SearchConstraints, SearchDiagnostics, SearchResult
    using .IsotopicPatternModule: IonSpec, ion_mass, IsotopePopulation, IsotopeContribution, IsotopePeak, PatternResult

    include("PatternMatching.jl")
    include("ResultExport.jl")
    export result_dict
    export PeakMatch, RankedCandidate, rank_formulas
    export IonSpec, ion_mass, IsotopePopulation, IsotopeContribution, IsotopePeak, PatternResult, SearchConstraints, SearchDiagnostics, SearchResult
    export isotopic_pattern, monoisotopic_mass, isotopic_pattern_protonated, monoisotopic_mass_protonated, find_formula, Compound

end # module
