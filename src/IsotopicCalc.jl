module IsotopicCalc

    using JSON

    const ELEMENTS = JSON.parsefile(Base.Filesystem.joinpath(dirname(@__FILE__),"elements.json"));

    include("functions.jl")

    export isotopicPattern, isotopicPatternProtonated, monoisotopicMassProtonated
    export isopat, isopatprot, monoisomass, monoisomassprot

end # module
