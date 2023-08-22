module IsotopicCalc

    using JSON

    include("functions.jl")

    export isotopicPattern, isotopicPatternProtonated, monoisotopicMassProtonated
    export isopat, isopatprot, monoisomass, monoisomassprot

end # module
