struct Compound
    formula::String
    adduct::String
    charge::Int
    mz::Float64
    ppm::Float64
end

const ADDUCTS = Dict("M+H" => 1.00782503223, "M+Na" => 22.989769282, "M+K" => 38.9637064864, "M-H" => -1.00782503223)

function generate_formulas(mz_input::Float64, atom_pool::Dict{String, Int}, adduct::String, charge::Int)
    formulas = []
    adduct_mass = get(ADDUCTS, adduct, 0.0)
    
    function generate_combinations(elements, counts, idx, current)
        if idx > length(elements)
            
            # Extract the counts of C, H, N, and O from 'current'
            C = current[findfirst(isequal("C"), elements)]
            H = current[findfirst(isequal("H"), elements)]
            N = current[findfirst(isequal("N"), elements)]
            O = current[findfirst(isequal("O"), elements)]
            
            # Apply heuristic rules
            if H >= 0.5C && O <= C
                M = sum(current[i] * ELEMENTS[elements[i]]["Relative Atomic Mass"][1] for i in eachindex(elements)) + adduct_mass - charge*0.0005485
                mz_calculated = charge == 0 ? M : M / abs(charge)
                formula = join([string(elements[i], current[i] == 1 ? "" : current[i]) for i in 1:length(elements) if current[i] > 0])
                ppm = ((mz_input - mz_calculated) / mz_calculated) * 1e6
                push!(formulas, Compound(formula, adduct, charge, mz_calculated, ppm))
            end
            return
        end
        for count in 0:counts[idx]
            current[idx] = count
            generate_combinations(elements, counts, idx + 1, current)
        end
    end
    
    elements = collect(keys(atom_pool))
    max_counts = []
    for elem in elements
        push!(max_counts,minimum([atom_pool[elem],round(Int, mz_input/ELEMENTS[elem]["Relative Atomic Mass"][1])]))
    end
    generate_combinations(elements, max_counts, 1, zeros(Int, length(elements)))
    
    return formulas
end

function filter_formulas(formulas, mz_input, tolerance)
    return filter(c -> abs(c.ppm) <= tolerance, formulas)
end

function findFormula(mz_input::Float64; tolerance_ppm::Number=100, atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct::String="", charge::Int=0)
    formulas = generate_formulas(mz_input, atom_pool, adduct, charge)
    matching_formulas = filter_formulas(formulas, mz_input, tolerance_ppm)
    
    # Sort the matching formulas by ppm
    sort!(matching_formulas, by = c -> abs(c.ppm))
    
    println("Matching formulas:")
    for compound in matching_formulas
        charge_sign = charge > 0 ? "+" : (charge < 0 ? "-" : "")
        adduct_string = adduct == "" ? charge_sign : "\t[" * adduct * "]" * charge_sign
        println(" ", compound.formula, adduct_string, "\tm/z: ", string(round(compound.mz, digits=6)), "\tppm: ", string(round(compound.ppm, digits=2)))
    end
    return matching_formulas
end