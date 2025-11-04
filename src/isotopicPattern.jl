# Pre-compiled regex patterns for better performance
const FORMULA_VALIDATION_REGEX = r"^[A-Za-z\[\]\d\(\)]+$"
const CHARGE_PATTERN_REGEX = r"\((\d+)([+-])\)"
const CHARGE_REMOVAL_REGEX = r"\(\d+[+-]\)"

# Parse the molecular formula
"""
    expand_round_brackets(formula::String) -> String

Expand parentheses in a chemical formula by multiplying the enclosed groups.

Recursively processes nested parentheses and applies multipliers.

# Arguments
- `formula::String`: Formula with parentheses, e.g., "(CH3)2CO" or "Ca(OH)2"

# Returns
- `String`: Expanded formula without parentheses, e.g., "CH3CH3CO" or "CaOH2"

# Examples
```julia
expand_round_brackets("(CH3)2CO") → "CH3CH3CO"
expand_round_brackets("Ca(OH)2") → "CaOHOH"
expand_round_brackets("((CH3)2)2") → "CH3CH3CH3CH3"
```
"""
function expand_round_brackets(formula::String)
    expanded_formula = ""
    i = 1
    n = length(formula)

    while i <= n
        c = formula[i]

        if c == '('
            j = i
            open_brackets = 1
            while j < n
                j += 1
                if formula[j] == '('
                    open_brackets += 1
                elseif formula[j] == ')'
                    open_brackets -= 1
                end
                if open_brackets == 0
                    break
                end
            end

            if open_brackets != 0
                throw(ArgumentError("Unmatched brackets in formula"))
            end

            subformula = formula[i+1:j-1]
            i = j + 1
            multiplier = ""

            while i <= n && isdigit(formula[i])
                multiplier *= formula[i]
                i += 1
            end

            multiplier = isempty(multiplier) ? 1 : parse(Int, multiplier)
            expanded_subformula = expand_round_brackets(subformula)  # recursively expand
            expanded_formula *= repeat(expanded_subformula, multiplier)
        else
            expanded_formula *= c
            i += 1
        end
    end
    return expanded_formula
end

function extract_square_brackets(formula::String)
    extracted = Dict{String, Int}()
    remaining_formula = ""
    i = 1
    len = length(formula)
    
    while i <= len
        c = formula[i]
        
        if c == '['  # Start of square bracket
            start_idx = i
            while i <= len && formula[i] != ']'
                i += 1
            end
            
            if i > len
                throw(ArgumentError("Mismatched square brackets in formula '$formula': opening '[' at position $start_idx has no closing ']'. Use square brackets for isotopes like [13C] or [2H]."))
            end

            element_inside_brackets = formula[start_idx+1:i-1]
            i += 1  # Move past the closing ']'
            
            # Check for a following multiplier
            multiplier = 1
            start_digit = i
            
            while i <= len && isdigit(formula[i])
                i += 1
            end
            
            if start_digit != i
                multiplier = parse(Int, formula[start_digit:i-1])
            end
            
            extracted[element_inside_brackets] = get(extracted, element_inside_brackets, 0) + multiplier
            
        else
            remaining_formula *= c
            i += 1
        end
    end
    
    return extracted, remaining_formula
end

function parse_formula_without_isotopes(remaining_formula::String)
    parsed = Dict{String, Int}()
    len = length(remaining_formula)
    i = 1
    current_element = ""
    while i <= len
        c = remaining_formula[i]
        if isuppercase(c)
            if current_element != ""
                parsed[current_element] = get(parsed, current_element, 0) + 1
            end
            current_element = string(c)
        elseif islowercase(c)
            current_element *= c
        elseif isdigit(c)
            start_idx = i
            while i <= len && isdigit(remaining_formula[i])
                i += 1
            end
            multiplier = parse(Int, remaining_formula[start_idx:i-1])
            parsed[current_element] = get(parsed, current_element, 0) + multiplier
            current_element = ""
            continue  # Skip the rest of the loop
        end
        i += 1
    end
    # If an element is left over at the end
    if current_element != ""
        parsed[current_element] = get(parsed, current_element, 0) + 1
    end
    return parsed
end

"""
    parse_formula(formula::String) -> Dict{String, Int}

Parse a chemical formula string into a dictionary of element counts.

Handles parentheses, square bracket isotope notation, and standard element notation.
This is the main formula parsing entry point.

# Arguments
- `formula::String`: Chemical formula (e.g., "C3H6O", "(CH3)2CO", "[13C]H3COCH3")

# Returns
- `Dict{String, Int}`: Dictionary mapping element/isotope symbols to their counts

# Examples
```julia
parse_formula("C3H6O") → Dict("C" => 3, "H" => 6, "O" => 1)
parse_formula("(CH3)2CO") → Dict("C" => 3, "H" => 6, "O" => 1)
parse_formula("[13C]H3COCH3") → Dict("13C" => 1, "C" => 2, "H" => 6, "O" => 1)
```
"""
function parse_formula(formula::String)
    expanded = expand_round_brackets(formula)
    extracted, remaining_formula = extract_square_brackets(expanded)
    if length(remaining_formula) > 0
        parsed_remaining = parse_formula_without_isotopes(remaining_formula)
        parsed_formula = merge(+, extracted, parsed_remaining)
    else
        parsed_formula = extracted
    end

    return parsed_formula
end

"""
    convolve(distro1, distro2, abundance_cutoff) -> Vector{Tuple{Float64, Float64}}

Convolve two isotopic distributions to combine them.

This implements the mathematical convolution operation for combining isotopic patterns.
Each peak in distro1 is combined with each peak in distro2 by adding masses and
multiplying abundances.

# Arguments
- `distro1`: First distribution as collection of (mass, abundance) tuples
- `distro2`: Second distribution as collection of (mass, abundance) tuples
- `abundance_cutoff`: Minimum abundance threshold to include in result

# Returns
- `Vector{Tuple{Float64, Float64}}`: Combined distribution sorted by mass

# Algorithm
For all combinations of peaks (m₁, a₁) and (m₂, a₂):
- Combined mass: m = m₁ + m₂
- Combined abundance: a = a₁ × a₂
- Include only if a ≥ abundance_cutoff
"""
function convolve(distro1, distro2, abundance_cutoff::Real)::Vector{Tuple{Float64, Float64}}
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

"""
    normalize_abundances(distribution) -> Vector{Tuple{Float64, Float64}}

Normalize isotopic abundances so the most abundant peak has abundance 1.0.

# Arguments
- `distribution`: Collection of (mass, abundance) tuples

# Returns
- `Vector{Tuple{Float64, Float64}}`: Normalized distribution where max abundance = 1.0
"""
function normalize_abundances(distribution::Vector{Tuple{Float64, Float64}})::Vector{Tuple{Float64, Float64}}
  # Find the maximum abundance
  max_abundance = maximum([a for (m, a) in distribution])

  # Normalize all abundances
  normalized_distribution = [(m, a / max_abundance) for (m, a) in distribution]

  return normalized_distribution
end

"""
    group_by_resolution(distribution, R) -> Vector{Tuple{Float64, Float64}}

Group isotopic peaks by mass resolution and sum their abundances.

Peaks within Δm/m = 1/R of each other are merged into a single peak.
This simulates the finite resolution of mass spectrometers.

# Arguments
- `distribution`: Collection of (mass, abundance) tuples
- `R`: Mass resolution (e.g., 10000 for R=10000)

# Returns
- `Vector{Tuple{Float64, Float64}}`: Distribution with peaks grouped by resolution

# Examples
At R=10000, peaks at m/z 100.0000 and 100.0050 (Δm=0.005) are separated
because Δm/m = 0.005/100 = 5×10⁻⁵ > 1/10000 = 1×10⁻⁴
"""
function group_by_resolution(distribution::Vector{Tuple{Float64, Float64}}, R::Real)::Vector{Tuple{Float64, Float64}}
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
  if occursin(CHARGE_PATTERN_REGEX, adduct)
      match_data = match(CHARGE_PATTERN_REGEX, adduct)
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

"""
    isotopicPattern(formula::String; abundance_cutoff=1e-5, R=10000, adduct::String="", print=true)

Calculate the isotopic pattern distribution for a given chemical formula.

This function computes all possible isotopic combinations and their relative abundances
using NIST isotopic composition data. The result includes both the mass and relative
abundance of each isotopologue.

# Arguments
- `formula::String`: Chemical formula in various formats:
  - Simple: "C3H6O"
  - With parentheses: "(CH3)2CO" or "CH3COCH3"
  - With specific isotopes: "[13C]H3COCH3" or "C3[2H]6O"
  - Deuterium shorthand: "C3D6O" (equivalent to C3[2H]6O)
- `abundance_cutoff::Number=1e-5`: Minimum relative abundance threshold (default: 10⁻⁵)
- `R::Number=10000`: Mass resolution for grouping nearby peaks (default: 10000)
- `adduct::String=""`: Optional adduct or ion (e.g., "H+", "Na+", "K+", "H-")
- `print::Bool=true`: Whether to print formatted output to console

# Returns
- `Vector{Tuple{Float64, Float64}}`: Array of (mass, abundance) tuples where:
  - mass: Isotopic mass in atomic mass units (amu)
  - abundance: Relative abundance (normalized to most abundant peak = 1.0)

# Examples
```julia
julia> isotopicPattern("CH3COCH3")
Formula:  C3H6O
----------------------------
Mass [amu]     Abundance [%]
----------------------------
58.0419        100.0000
59.0452        3.2447
...

julia> isotopicPattern("C6H12O6"; abundance_cutoff=1e-3, adduct="Na+")
# Returns distribution for sodium adduct of glucose with 0.1% cutoff

julia> pattern = isotopicPattern("C3H6O"; print=false)
# Returns data without printing (useful for programmatic access)
```

# Notes
- Uses convolution algorithm to combine isotopic distributions
- Resolution parameter `R` groups peaks within Δm/m = 1/R
- Electron mass adjustments applied for charged species
- Based on NIST Atomic Weights and Isotopic Compositions database

# Throws
- `ArgumentError`: If formula contains invalid characters, brackets are unbalanced,
  or if abundance_cutoff < 0 or R ≤ 0
- `ArgumentError`: If unknown element or isotope is specified

See also: [`monoisotopicMass`](@ref), [`isotopicPatternProtonated`](@ref), [`findFormula`](@ref)
"""
function isotopicPattern(formula::String; abundance_cutoff=1e-5, R=10000, adduct::String="", print=true)  # Default resolution R = 10000
  # Validate formula format, abundance cutoff, and mass resolution
  if isempty(strip(formula))
    throw(ArgumentError("Formula cannot be empty. Provide a chemical formula like 'H2O', 'C6H12O6', or 'CH3COOH'."))
  end
  if !occursin(FORMULA_VALIDATION_REGEX, formula)
    throw(ArgumentError("Invalid formula '$formula'. Use only letters, numbers, parentheses (), and square brackets []. Examples: 'C3H6O', '(CH3)2CO', '[13C]H3COCH3'"))
  end
  if abundance_cutoff < 0
    throw(ArgumentError("abundance_cutoff must be non-negative (got $abundance_cutoff). Use values like 1e-5 for 0.001% cutoff."))
  end
  if R <= 0
    throw(ArgumentError("Resolution R must be positive (got $R). Typical values: 1000-100000 (e.g., R=10000 for FT-ICR)."))
  end
  if !(count(x -> x == '(', formula) == count(x -> x == ')', formula))
    open_count = count(x -> x == '(', formula)
    close_count = count(x -> x == ')', formula)
    throw(ArgumentError("Unbalanced parentheses in formula '$formula': $open_count opening '(' but $close_count closing ')'. Each '(' must have a matching ')'."))
  end
  parsed_formula = parse_formula(formula)
  # Cache the formatted formula string to avoid duplicate parsing
  formula_display = sum_formula(parsed_formula)

  # If an adduct is provided, parse and combine it with the base formula
  if !isempty(adduct)
      parsed_adduct = parse_formula(replace(adduct, CHARGE_REMOVAL_REGEX => ""))  # Remove charge info for parsing
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
              # Suggest common elements
              common_elements = ["H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "Na", "K", "Ca", "Fe"]
              throw(ArgumentError("Unknown element '$base_element' (from '$atom'). Common elements: $(join(common_elements, ", ")). Check spelling and capitalization."))
          end
          isotope_index = findfirst(x -> x == atom, ELEMENTS[base_element]["Symbol"])
          if isnothing(isotope_index)
              available_isotopes = join(ELEMENTS[base_element]["Symbol"], ", ")
              throw(ArgumentError("Unknown isotope '$atom'. Available isotopes of $base_element: $available_isotopes. Use bracket notation like [13C] or [2H]."))
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
    output = "\nFormula: \e[1;32m $formula_display"
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
        m_round = Printf.Format("%."*string(maximum([3, round(Int, 1-log10((m/R)/3))]))*"f")
        mIO = IOBuffer(); Printf.format(mIO, m_round, m); mm = String(take!(mIO))

        # abbundance in % formated to 3 decimal point even padded with zero if necessary
        a_round = Printf.Format("%."*string(round(Int,1+log(10,1/abundance_cutoff/100)))*"f")
        aIO = IOBuffer(); Printf.format(aIO, a_round, a*100); aa = String(take!(aIO))

        # padding list dependent on mass
        format_pad = Printf.Format("%"*string(15-length(mm)+length(aa))*"s")
        padIO = IOBuffer(); Printf.format(padIO, format_pad, aa); pad = String(take!(padIO))

        output *= mm * pad *"\n"
    end
    output *= "Found $(length(final_distribution)) isotopic masses for $abundance_cutoff abundance limit."

    println(output)
  end
  return final_distribution
end

# special cases
"""
    isotopicPatternProtonated(formula::String)

Convenience function to calculate the isotopic pattern for a protonated molecule [M+H]⁺.

This is equivalent to calling `isotopicPattern(formula; adduct="H+")`.

# Arguments
- `formula::String`: Chemical formula (same format as `isotopicPattern`)

# Returns
- `Vector{Tuple{Float64, Float64}}`: Array of (mass, abundance) tuples for [M+H]⁺

# Examples
```julia
julia> isotopicPatternProtonated("C3H6O")
# Calculates pattern for protonated acetone (m/z 59.049)
```

See also: [`isotopicPattern`](@ref), [`monoisotopicMassProtonated`](@ref)
"""
isotopicPatternProtonated(x) = isotopicPattern(x; adduct="H+")

"""
    monoisotopicMass(formula::String)

Calculate the monoisotopic mass (most abundant isotopic composition) for a given formula.

Returns only the mass of the lowest-mass isotopologue (all lightest isotopes),
without calculating the full isotopic distribution or printing output.

# Arguments
- `formula::String`: Chemical formula (same format as `isotopicPattern`)

# Returns
- `Float64`: Monoisotopic mass in atomic mass units (amu)

# Examples
```julia
julia> monoisotopicMass("CH3COCH3")
58.0419

julia> monoisotopicMass("C6H12O6")
180.0634

julia> monoisotopicMass("[13C]H3COCH3")
59.0452  # One carbon is ¹³C instead of ¹²C
```

# Notes
- This is computationally efficient as it only extracts the first mass peak
- Uses natural isotopic abundances from NIST database
- For charged species, use the `adduct` parameter in `isotopicPattern` instead

See also: [`isotopicPattern`](@ref), [`monoisotopicMassProtonated`](@ref)
"""
monoisotopicMass(x) = isotopicPattern(x;print=false)[1][1]

"""
    monoisotopicMassProtonated(formula::String)

Calculate the monoisotopic mass of a protonated molecule [M+H]⁺.

Equivalent to calling `isotopicPattern(formula; adduct="H+", print=false)[1][1]`.
Useful for mass spectrometry applications where positive-mode ESI is common.

# Arguments
- `formula::String`: Chemical formula (same format as `isotopicPattern`)

# Returns
- `Float64`: Monoisotopic m/z value for [M+H]⁺ in atomic mass units

# Examples
```julia
julia> monoisotopicMassProtonated("C3H6O")
59.049141  # Acetone + proton

julia> monoisotopicMassProtonated("C6H12O6")
181.070646  # Glucose + proton
```

See also: [`monoisotopicMass`](@ref), [`isotopicPatternProtonated`](@ref)
"""
monoisotopicMassProtonated(x) = isotopicPattern(x; adduct="H+", print=false)[1][1]
