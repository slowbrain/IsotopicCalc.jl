module IsotopicCalc

    using JSON, Printf

    const ELEMENTS = JSON.parsefile(Base.Filesystem.joinpath(dirname(@__FILE__),"elements.json"));

    include("isotopicPattern.jl")
    include("findFormula.jl")

    export isotopicPattern, monoisotopicMass, isotopicPatternProtonated, monoisotopicMassProtonated, findFormula, Compound

end # module
