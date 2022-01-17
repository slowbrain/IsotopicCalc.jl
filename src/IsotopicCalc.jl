module IsotopicCalc

    using JSON, OrderedCollections

    include("functions.jl")

    export isotopicPattern, isotopicPatternProtonated, monoisotopicMass
    export isopat, isopatprot, monoisomass

end # module
