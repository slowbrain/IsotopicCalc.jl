struct Compound
    formula::String
    adduct::String
    charge::Int
    mz::Float64
    ppm::Float64
end

const ATOMIC_MASSES = Dict("C" => 12.0, "H" => 1.00782503223, "N" => 14.00307400443, "O" => 15.99491461957)
const ADDUCTS = Dict("M+H" => 1.00782503223, "M+Na" => 22.989769282, "M+K" => 38.9637064864, "M-H" => -1.00782503223)

function generate_compounds(atom_pool::Dict{String, Int}, adduct::String, charge::Int, mz_input::Float64)
    compounds = []
    adduct_mass = get(ADDUCTS, adduct, 0.0)
    
    function generate_combinations(elements, counts, idx, current)
        if idx > length(elements)
            
            # Extract the counts of C, H, N, and O from 'current'
            C = current[findfirst(isequal("C"), elements)]
            H = current[findfirst(isequal("H"), elements)]
            N = current[findfirst(isequal("N"), elements)]
            O = current[findfirst(isequal("O"), elements)]
            
            # Apply heuristic rules
            if N <= 0.5C
                M = sum(current[i] * ATOMIC_MASSES[elements[i]] for i in eachindex(elements)) + adduct_mass - charge*0.0005485
                mz_calculated = charge == 0 ? M : M / abs(charge)
                formula = join([string(elements[i], current[i] == 1 ? "" : current[i]) for i in 1:length(elements) if current[i] > 0])
                ppm = ((mz_input - mz_calculated) / mz_calculated) * 1e6
                push!(compounds, Compound(formula, adduct, charge, mz_calculated, ppm))
            end
            return
        end
        for count in 0:counts[idx]
            current[idx] = count
            generate_combinations(elements, counts, idx + 1, current)
        end
    end
    
    elements = collect(keys(atom_pool))
    max_counts = collect(values(atom_pool))
    generate_combinations(elements, max_counts, 1, zeros(Int, length(elements)))
    
    return compounds
end

function filter_compounds(compounds, mz_input, tolerance)
    return filter(c -> abs(c.ppm) <= tolerance, compounds)
end

function find_compounds(mz_input::Float64, tolerance::Float64, atom_pool::Dict{String, Int}, adduct::String, charge::Int)
    compounds = generate_compounds(atom_pool, adduct, charge, mz_input)
    matching_compounds = filter_compounds(compounds, mz_input, tolerance)
    
    println("Matching Compounds:")
    for compound in matching_compounds
        charge_sign = charge > 0 ? "+" : (charge < 0 ? "-" : "")
        # adduct_element = replace(adduct, "M" => "")
        # adduct_string = adduct == "" ? charge_sign : "." * adduct_element * charge_sign
        adduct_string = adduct == "" ? charge_sign : "\t[" * adduct * "]" * charge_sign
        println(" ", compound.formula, adduct_string, "\tm/z: ", string(round(compound.mz, digits=6)), "\tppm: ", string(round(compound.ppm, digits=2)))
    end
    return matching_compounds
end