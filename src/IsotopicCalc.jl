module IsotopicCalc

    using JSON, Printf

    const ELEMENTS = JSON.parsefile(Base.Filesystem.joinpath(dirname(@__FILE__),"elements.json"));

    include("isotopicPattern.jl")
    include("determineFormula.jl")

    export isotopicPattern, isotopicPatternProtonated, monoisotopicMassProtonated
    export isopat, isopatprot, monoisomass, monoisomassprot

end # module
