function cart(x::Array,y::Array)
    #return [repmat(x,1,length(y))'[:] repmat(y,length(x),1)[:]]
    leny=length(y)
    lenx=length(x)
    m=leny*lenx
    OUT = zeros(Float64, m, 2)
    c=1
    for i = 1:lenx
        for j = 1:leny
            OUT[c,1] = x[i]
            OUT[c,2] = y[j]
            c+=1
        end
    end
    return OUT
end

function cartSum(x::Array,y::Array)
    return sum(cart(x, y), dims=2)
end

function cartProd(x::Array,y::Array)
    return prod(cart(x, y), dims=2)
end

function parseFormula(formula::AbstractString)
  # support for any kind of brackets
  formula = replace(formula, "[" => "(")
  formula = replace(formula, "]" => ")")
  formula = replace(formula, "{" => "(")
  formula = replace(formula, "}" => ")")
  # check for any existing brackets
  formula_no_brackets = formula
  while occursin(r"(.*)\((.*?)\)(\d*)", formula_no_brackets)
    # find rightmost left parenthesis
    formula_split = match(r"(.*)\((.*?)\)(\d*)", formula_no_brackets);
    leftpart = formula_split.captures[1];
    bracket_content = formula_split.captures[2];
    if (formula_split.captures[3] == "")
      bracket_content_mult = 1
    else
      bracket_content_mult = parse(Int,formula_split.captures[3])
    end
    rightpart = replace(formula_no_brackets, formula_split.match => "")
    formula_no_brackets = leftpart*repeat(bracket_content, bracket_content_mult)*rightpart
  end

  elementsRawIDs = collect((m.match for m = eachmatch(r"([A-Z][a-z]{0,1})(\d*)", string(formula_no_brackets))))
  #elementsRawIDs = matchall(r"([A-Z](?:[a-z])?)(\d*)", string(formula))

  elementsDict = Dict{String, Int}()
  for i in elementsRawIDs
    if occursin(r"[0-9]{1,}", string(i))
      elementID = match(r"[A-Za-z]*", i).match
      elementMult = parse(Int, match(r"\d*$", i).match)
      try
        elementsDict[elementID] = elementsDict[elementID] + elementMult
      catch
        elementsDict[elementID] = elementMult
      end
    else
      try
        elementsDict[i] = elementsDict[i] + 1
      catch
        elementsDict[i] = 1
      end
    end
  end
  return elementsDict
end

function parseCharge(charge::AbstractString)
  # support for any kind of brackets
  charge = replace(charge, "[" => "(")
  charge = replace(charge, "]" => ")")
  charge = replace(charge, "{" => "(")
  charge = replace(charge, "}" => ")")

  if occursin(r"(\d)\((\+*\-*)\)", charge)
    # find rightmost left parenthesis
    charge_split = match(r"(\d)\((\+*\-*)\)", charge);
    charge_multiple = charge_split.captures[1];
    charge_polarity = charge_split.captures[2];
  else
    charge_split = match(r"(\d)*(\+*\-*)", charge);
    charge_multiple = charge_split.captures[1];
    charge_polarity = charge_split.captures[2];
  end

  return (charge_multiple, charge_polarity)
end


function convoluteDict(sumForm::Dict)
  cD = Array{String,1}()
  for i in keys(sumForm)
    for j in 1:sumForm[i]
      push!(cD, i)
    end
  end
  return cD
end

function elementsStringify!(stringBase::AbstractString, a::AbstractString, b::Number)
  if b == 1
    stringBase = stringBase*string(a)
  else
    stringBase = stringBase*string(a)*string(b)
  end
  stringBase
end

function sumFormula(sumFormDict::Dict)
  sf = ""
  sf_c_prefix = ""
  sf_h_prefix = ""
  for i in keys(sumFormDict)
    if string(i) == "C"
      sf_c_prefix = elementsStringify!(sf_c_prefix, "C", sumFormDict[i])
    elseif string(i) == "H"
      sf_h_prefix = elementsStringify!(sf_h_prefix, "H", sumFormDict[i])
    else
      sf = elementsStringify!(sf, string(i), sumFormDict[i])
    end
  end
  return sf_c_prefix*sf_h_prefix*sf
end

"""
Usage: isotopicPattern(formula::AbstractString; adduct=\"H\", charge=\"+\", abundanceLimit = 1e-5, niceOutput = true)")

args:   formula - can be sum formula
kwargs: adduct - possible adduct to the  
"""
function isotopicPattern(formula::String; kwargs...)
    kwargs = Dict(kwargs)
    if haskey(kwargs, :adduct)
      adduct = kwargs[:adduct]
    else
      adduct = ""
    end

    if haskey(kwargs, :charge)
      charge = kwargs[:charge]
    else
      charge = ""
    end

    if haskey(kwargs, :abundanceLimit)
      abundanceLimit = kwargs[:abundanceLimit]
    else
      abundanceLimit = 1e-5
    end

    if haskey(kwargs, :niceOutput)
      niceOutput = kwargs[:niceOutput]
    else
      niceOutput = true
    end

    if haskey(kwargs, :result)
      exportBool = kwargs[:result]
      if exportBool
        niceOutput = false
      end
    else
      exportBool = false
    end

    start =  Base.time()
    elementsDict = parseFormula(formula)
    composition = convoluteDict(elementsDict)
    #println(composition)

    if length(adduct) > 0
        composition = [composition; convoluteDict(parseFormula(adduct))]
    end

    elem = string(composition[1])
    masses = ELEMENTS[elem]["Relative Atomic Mass"]
    abund = ELEMENTS[elem]["Isotopic Composition"]

    for i in 2:length(composition)
        elemNext = string(composition[i])
        massesNext = ELEMENTS[elemNext]["Relative Atomic Mass"]
        abundNext = ELEMENTS[elemNext]["Isotopic Composition"]

        massesTemp = round.(cartSum(masses, massesNext), digits = 6)
        abundTemp = cartProd(abund, abundNext)

        masses = unique(massesTemp)

        abund = Array{Float64}(undef, size(masses))
        for j in 1:length(masses)
            index = findall((in)(masses[j]), massesTemp)
            abund[j] = sum(abundTemp[index])
        end

        toKeep = findall(abund.>abundanceLimit)
        masses = masses[toKeep]
        abund = abund[toKeep]
    end
    #println(abund)

    SF = sumFormula(elementsDict)
    SF_add = ""
    if length(adduct) > 0
      SF_add = ".$(sumFormula(parseFormula(adduct)))"
    end

    if length(charge) > 0
      (charge_multiple, charge_polarity) = parseCharge(charge)

      if charge_polarity in ["+", "-"]
        try
          parse(Int,charge_multiple)
        catch
          charge_multiple = "1"
        end
        for i in 1:parse(Int,charge_multiple)
          masses = [eval(Meta.parse(charge_polarity))(masses[j], -1*ELEMENTS["e"]["Relative Atomic Mass"][1]) for j in 1:length(masses)]
        end
        if parse(Int,charge_multiple) > 1
          masses = masses/parse(Int,charge_multiple)
        end
      end
    else
      charge_multiple = ""
      charge_polarity = ""
    end

    charge_multiple_label = ""
    if charge_multiple != "1"
      try
        charge_multiple_label = Base.REPLCompletions.latex_symbols["\\^$(charge_multiple)"]
      catch
        charge_multiple_label = charge_multiple
      end
    else
      charge_multiple_label = ""
    end

    charge_polarity_label = ""
    try
      charge_polarity_label = Base.REPLCompletions.latex_symbols["\\^$(charge_polarity)"]
    catch
      charge_polarity_label = charge_polarity
    end

    SFx = SF*SF_add*charge_multiple_label*charge_polarity_label

    if niceOutput
      println("\nFormula: \e[1;32m $SFx")
    end

    IPCresult = Dict( SFx => zeros(Float64, length(masses), 2))
    #println(result)
    IPCresult[SFx][:,2] = abund/maximum(abund)
    IPCresult[SFx][:,1] = masses

    IPCresult[SFx] = IPCresult[SFx][sortperm(IPCresult[SFx][:,2], rev=true),:]

    IPCresult[SFx] = IPCresult[SFx][findall(IPCresult[SFx][:,2].>abundanceLimit),:]
    stop = time()
    #println("\e[1;37mevaluated in ", round(stop-start,3), " seconds." )
    if niceOutput
      println("\e[m-------------------------------")
      println("Mass [amu]\tAbbundance [%]")
      println("-------------------------------")
      for i in 1:size(IPCresult[SFx])[1]
          println(string(round(IPCresult[SFx][i,1], digits=6))*"\t"*string(round(IPCresult[SFx][i,2]*100,digits=6)))
        end
        println("Found ", length(masses), " isotopic masses for $abundanceLimit abundance limit.")
        println()
      end
    # if exportBool
    #   return IPCresult[SFx]
    # end
    return IPCresult[SFx]
end;

isopat(x...) = isotopicPattern(x...)

# special cases
isotopicPatternProtonated(x) = isotopicPattern(x; adduct="H", charge="+")
isopatprot(x...) = isotopicPatternProtonated(x...)

monoisotopicMass(x) = isotopicPattern(x...; niceOutput=false)[1,1]
monoisomass(x...) = monoisotopicMass(x...)

monoisotopicMassProtonated(x) = isotopicPattern(x...; adduct="H", charge="+", niceOutput=false)[1,1]
monoisomassprot(x...) = monoisotopicMassProtonated(x...)



# Parse the molecular formula
function parse_formula(formula::String)
  atoms = Dict{String, Int}()
  
  # Handle specific isotopes in square brackets
  isotopic_matches = eachmatch(r"\[([^\]]+)\]", formula)
  for m in isotopic_matches
      isotope = m.captures[1]
      atoms[isotope] = get(atoms, isotope, 0) + 1
  end
  formula = replace(formula, r"\[([^\]]+)\]" => "")
  
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
      rounded_mass = round(m, digits=maximum(3, round(Int, -log10(Δm/10))))
      
      # Group by the rounded mass and sum the abundances
      grouped[rounded_mass] = get(grouped, rounded_mass, 0.0) + a
  end
  
  return sort(collect(grouped))
end

function isotopic_pattern(formula::String, abundance_cutoff=1e-5, R=5000)  # Default resolution R = 5000
  parsed_formula = parse_formula(formula)
  
  # Initial distribution: monoisotopic mass with 100% abundance
  final_distribution = [(0.0, 1.0)]
  
  for (atom, count) in parsed_formula
      masses = ELEMENTS[atom]["Relative Atomic Mass"]
      abundances = ELEMENTS[atom]["Isotopic Composition"]
      atom_distribution = sort(collect(zip(masses, abundances)))
      
      for _ = 1:count
          final_distribution = convolve(final_distribution, atom_distribution, abundance_cutoff)
          # Filter out isotopes below the abundance cutoff
          final_distribution = filter(x -> x[2] >= abundance_cutoff, final_distribution)
      end
  end
  
  # Group and sum abundances based on the provided mass resolution R
  final_distribution = group_by_resolution(final_distribution, R)
  
  # Normalize the abundances to 100%
  final_distribution = normalize_abundances(final_distribution)
  
  # Sort the distribution by mass
  return sort(final_distribution)
end
