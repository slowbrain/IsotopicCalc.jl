# Parse the molecular formula
function parse_formula(formula::String)
  parsed_formula = Dict{String, Int}()
  i = 1
  n = length(formula)
  local_multiplier = 1

  while i <= n
      c = formula[i]
      if isletter(c)
          k = i + 1
          while k <= n && (isletter(formula[k]) || isdigit(formula[k]))
              k += 1
          end
          atom = formula[i:k-1]
          count = 1  # Default count is 1
          if k <= n && isdigit(formula[k])
              l = k
              while l <= n && isdigit(formula[l])
                  l += 1
              end
              count = parse(Int, formula[k:l-1])
              k = l
          end
          parsed_formula[atom] = get(parsed_formula, atom, 0) + (count * local_multiplier)
          i = k
      elseif c == '(' || c == '['
          closing = (c == '(' ? ')' : ']')
          counter = 1
          k = i + 1
          while k <= n
              if formula[k] == c
                  counter += 1
              elseif formula[k] == closing
                  counter -= 1
              end
              if counter == 0
                  break
              end
              k += 1
          end
          if counter != 0
              throw(ArgumentError("Unmatched $c in formula"))
          end
          sub_formula = formula[i+1:k-1]
          k += 1  # Move past the closing bracket
          multiplier = 1  # Default multiplier is 1
          if k <= n && isdigit(formula[k])
              l = k
              while l <= n && isdigit(formula[l])
                  l += 1
              end
              multiplier = parse(Int, formula[k:l-1])
              k = l
          end
          sub_parsed = parse_formula(sub_formula)
          for (atom, count) in sub_parsed
              parsed_formula[atom] = get(parsed_formula, atom, 0) + (count * multiplier * local_multiplier)
          end
          i = k
      elseif isdigit(c)
          i += 1
      else
          throw(ArgumentError("Invalid character $c in formula"))
      end
  end
  return parsed_formula
end

parsed_formula = parse_formula("([13C]H3)2CO")
println(parsed_formula)

function convolve(distro1, distro2, abundance_cutoff)
  new_distribution = Dict{Float64, Float64}()
  sorted1 = sort(collect(distro1))
  sorted2 = sort(collect(distro2))
  
  for (m1, a1) in sorted1
      for (m2, a2) in sorted2
          mass = m1 + m2
          abundance = a1 * a2
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
      rounded_mass = round(m, digits=maximum([3, round(Int, 1-log10(Δm/3))]))
      
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

function sum_formula(parsed_formula::Dict{String, Int})
  elements_order = ["C", "H", "O", "N", "P", "S", "F", "Cl", "Br", "I"]  # Common elements first
  other_elements = setdiff(collect(keys(parsed_formula)), elements_order)  # Convert keys to an array
  sort!(other_elements)

  sum_formula = ""
  for element in [elements_order; other_elements]
      count = get(parsed_formula, element, 0)
      if count > 0
          sum_formula *= element
          sum_formula *= count > 1 ? string(count) : ""
      end
  end
  return sum_formula
end


function combine_formulas(base_formula::Dict{String, Int}, adduct_formula::Dict{String, Int})
  combined_formula = deepcopy(base_formula)
  for (atom, count) in adduct_formula
      combined_formula[atom] = get(combined_formula, atom, 0) + count
  end
  return combined_formula
end

function extract_charge(adduct::String)
  # Special case for single "+" or "-"
  if adduct == "+" return -1 end
  if adduct == "-" return 1 end

  # Handle other formats with regex
  if occursin(r"\((\d+)([+-])\)", adduct)
      match_data = match(r"\((\d+)([+-])\)", adduct)
      charge_num = parse(Int, match_data.captures[1])
      charge_sign = match_data.captures[2]
      return charge_sign == "+" ? -charge_num : charge_num  # Return negative for cationic and positive for anionic
  elseif endswith(adduct, "+")
      return -1  # cationic
  elseif endswith(adduct, "-")
      return 1   # anionic
  else
      return 0  # No charge information found
  end
end

function isotopicPattern(formula::String; abundance_cutoff=1e-5, R=10000, adduct::String="", print=true)  # Default resolution R = 10000
  # Validate formula format, abundance cutoff, and mass resolution
  if !occursin(r"^([A-Z][a-z]*\d*)*(\([\w]+\)\d*)*$", formula) || abundance_cutoff < 0 || R <= 0
    throw(ArgumentError("Invalid input parameters"))
  end
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
    if atom != "e"  # Skip electrons
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

  if print
    # Format the output
    output = "\nFormula: \e[1;32m $(sum_formula(parse_formula(formula)))"
    if !isempty(adduct)
        output *= ".$adduct"
    end
    output *= "\n\e[m----------------------------\n"
    # constant padding in title independent on tab length
    pad_tileIO = IOBuffer(); @printf(pad_tileIO, "%5s", ""); pad_title = String(take!(pad_tileIO))
    output *= "Mass [amu]"*pad_title*"Abundance [%]\n"
    output *= "----------------------------\n"
    for (m, a) in final_distribution
        # mass formated to 4 decimal point even padded with zero if necessary
        # mIO = IOBuffer(); @printf(mIO, "%.4f", m); mm = String(take!(mIO))
        m_round = Printf.Format("%."*string(maximum([3, round(Int, 1-log10((m/R)/3))]))*"f")
        mIO = IOBuffer(); Printf.format(mIO, m_round, m); mm = String(take!(mIO))

        # abbundance in % formated to 3 decimal point even padded with zero if necessary
        # aIO = IOBuffer(); @printf(aIO, "%.3f", a*100); aa = String(take!(aIO))
        a_round = Printf.Format("%."*string(round(Int,1+log(10,1/abundance_cutoff/100)))*"f")
        aIO = IOBuffer(); Printf.format(aIO, a_round, a*100); aa = String(take!(aIO))

        # padding list dependent on mass
        format_pad = Printf.Format("%"*string(15-length(mm)+length(aa))*"s")
        padIO = IOBuffer(); Printf.format(padIO, format_pad, aa); pad = String(take!(padIO))
        #padIO = IOBuffer(); @printf(padIO, "%14s", aa); pad = String(take!(padIO))

        output *= mm * pad *"\n"
    end
    output *= "Found $(length(final_distribution)) isotopic masses for $abundance_cutoff abundance limit."

    println(output)
  end
  return final_distribution
end

# isopat(x...) = isotopicPattern(x...)

# special cases
isotopicPatternProtonated(x) = isotopicPattern(x; adduct="H+")

monoisotopicMass(x) = isotopicPattern(x;print=false)[1][1]

monoisotopicMassProtonated(x) = isotopicPattern(x; adduct="H+", print=false)[1][1]
