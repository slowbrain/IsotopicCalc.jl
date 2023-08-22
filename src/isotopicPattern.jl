# Parse the molecular formula
function parse_formula(formula::String)
  atoms = Dict{String, Int}()
  
  # Handle specific isotopes in square brackets
  isotopic_matches = eachmatch(r"\[([^\]]+)\](\d*)", formula)
  for m in isotopic_matches
      isotope = m.captures[1]
      count_str = m.captures[2]
      count = isempty(count_str) ? 1 : parse(Int, count_str)
      atoms[isotope] = get(atoms, isotope, 0) + count
  end
  formula = replace(formula, r"\[([^\]]+)\](\d*)" => "")
  
  # Expand bracketed groups
  while true
      match_data = match(r"\(([^)]+)\)(\d+)", formula)
      if match_data === nothing
          break
      end
      group_content = match_data.captures[1]
      multiplier = parse(Int, match_data.captures[2])
      replacement = repeat(group_content, multiplier)
      formula = replace(formula, match_data.match => replacement)
  end
  
  # Match element symbols followed by optional digit counts
  matches = eachmatch(r"([A-Z][a-z]*)(\d*)", formula)
  
  for m in matches
      element = m.captures[1]
      count_str = m.captures[2]
      count = isempty(count_str) ? 1 : parse(Int, count_str)
      
      # Add the count of the element to the dictionary, accumulating if the element already exists
      atoms[element] = get(atoms, element, 0) + count
  end
  
  return atoms
end

function convolve(distro1, distro2, abundance_cutoff)
  new_distribution = Dict{Float64, Float64}()
  for (m1, a1) in distro1
      for (m2, a2) in distro2
          mass = m1 + m2
          abundance = a1 * a2
          # Only include isotopes with abundances above the cutoff
          if abundance >= abundance_cutoff
              new_distribution[mass] = get(new_distribution, mass, 0.0) + abundance
          end
      end
  end
  return sort(collect(new_distribution))
end

function normalize_abundances(distribution)
  # Find the maximum abundance
  max_abundance = maximum([a for (m, a) in distribution])
  
  # Normalize all abundances
  normalized_distribution = [(m, a / max_abundance) for (m, a) in distribution]
  
  return normalized_distribution
end

function group_by_resolution(distribution, R)
  grouped = Dict{Float64, Float64}()
  
  for (m, a) in distribution
      Δm = m / R
      # Round the mass based on Δm and ensure a minimum of 3 decimal places
      rounded_mass = round(m, digits=maximum([4, round(Int, -log10(Δm/10))]))
      
      # Group by the rounded mass and sum the abundances
      grouped[rounded_mass] = get(grouped, rounded_mass, 0.0) + a
  end
  
  return sort(collect(grouped))
end

function get_base_element(isotope::String)
  if isuppercase(isotope[end])
      return string(isotope[end])
  else
      return isotope[end-1:end]
  end
end

function combine_formulas(base_formula::Dict{String, Int}, adduct_formula::Dict{String, Int})
  combined_formula = deepcopy(base_formula)
  for (atom, count) in adduct_formula
      combined_formula[atom] = get(combined_formula, atom, 0) + count
  end
  return combined_formula
end

function extract_charge(adduct::String)
  # Handle format "Ca(2+)" or "Cl(1-)"
  match_data = match(r"\((\d+)([+-])\)", adduct)
  if match_data !== nothing
      charge_num = parse(Int, match_data.captures[1])
      charge_sign = match_data.captures[2]
      return charge_sign == "+" ? -charge_num : charge_num  # Return negative for cationic and positive for anionic
  end

  # Handle format "H+" or "Na+"
  if endswith(adduct, "+")
      return -1  # cationic
  elseif endswith(adduct, "-")
      return 1   # anionic
  end

  return 0  # No charge information found
end

function isotopicPattern(formula::String; abundance_cutoff=1e-5, R=4000, adduct::String="")  # Default resolution R = 10000
  parsed_formula = parse_formula(formula)
    
  # If an adduct is provided, parse and combine it with the base formula
  if !isempty(adduct)
      parsed_adduct = parse_formula(replace(adduct, r"\(\d+[+-]\)" => ""))  # Remove charge info for parsing
      parsed_formula = combine_formulas(parsed_formula, parsed_adduct)
      
      # Extract charge from adduct and adjust electron count
      charge = extract_charge(adduct)
      parsed_formula["e"] = get(parsed_formula, "e", 0) + charge
  end

  # Initial distribution: monoisotopic mass with 100% abundance
  final_distribution = [(0.0, 1.0)]
  
  for (atom, count) in parsed_formula
      if haskey(ELEMENTS, atom)
          masses = ELEMENTS[atom]["Relative Atomic Mass"]
          abundances = ELEMENTS[atom]["Isotopic Composition"]
      else
          # Handle specific isotopes
          base_element = get_base_element(atom)
          if !haskey(ELEMENTS, base_element)
              throw(ArgumentError("Unknown element $base_element derived from isotope $atom"))
          end
          isotope_index = findfirst(x -> x == atom, ELEMENTS[base_element]["Symbol"])
          if isnothing(isotope_index)
              throw(ArgumentError("Unknown isotope $atom found in formula"))
          end
          masses = [ELEMENTS[base_element]["Relative Atomic Mass"][isotope_index]]
          abundances = [1.0]  # Specific isotopes are considered to have 100% abundance
      end
      
      atom_distribution = sort(collect(zip(masses, abundances)))
      
      for _ = 1:count
          final_distribution = convolve(final_distribution, atom_distribution, abundance_cutoff)
          # Filter out isotopes below the abundance cutoff
          final_distribution = filter(x -> x[2] >= abundance_cutoff, final_distribution)
      end
  end
  
  electron_count = get(parsed_formula, "e", 0)
  electron_mass = ELEMENTS["e"]["Relative Atomic Mass"][1]
  final_distribution = [(m + electron_count * electron_mass, a) for (m, a) in final_distribution]

  # Group and sum abundances based on the provided mass resolution R
  final_distribution = group_by_resolution(final_distribution, R)
  
  # Normalize the abundances to 100%
  final_distribution = normalize_abundances(final_distribution)
  
  # Sort the distribution by mass
  sort!(final_distribution)

  # Format the output
  output = "\nFormula: \e[1;32m $formula"
  if !isempty(adduct)
      output *= ".$adduct"
  end
  output *= "\n\e[m-------------------------------\n"
  output *= "Mass [amu]\tAbundance [%]\n"
  output *= "-------------------------------\n"
  for (m, a) in final_distribution
      if a >= 1
        output *= "$m\t\t$(round(Int,a * 100)).000\n"
      elseif a >= 0.1 && a < 1
        output *= "$m\t\t $(round(a * 100, digits=3))\n"
      elseif a < 0.1
        output *= "$m\t\t  $(round(a * 100, digits=3))\n"
      end
  end
  output *= "Found $(length(final_distribution)) isotopic masses for $abundance_cutoff abundance limit."
  
  println(output)
  return final_distribution
end

isopat(x...) = isotopicPattern(x...)

# special cases
isotopicPatternProtonated(x) = isotopicPattern(x; adduct="H+")
isopatprot(x...) = isotopicPatternProtonated(x...)

monoisotopicMass(x) = isotopicPattern(x...)[1][1]
monoisomass(x...) = monoisotopicMass(x...)

monoisotopicMassProtonated(x) = isotopicPattern(x...; adduct="H+")[1][1]
monoisomassprot(x...) = monoisotopicMassProtonated(x...)