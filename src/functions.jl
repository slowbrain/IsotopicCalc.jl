elementsAbundances = JSON.parsefile(Base.Filesystem.joinpath(dirname(@__FILE__),"elements.json"))
    
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
  return OrderedDict(elementsDict)
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

function convoluteDict(sumForm::OrderedDict)
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

function sumFormula(sumFormDict::OrderedDict)
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
    masses = elementsAbundances[elem]["Relative Atomic Mass"]
    abund = elementsAbundances[elem]["Isotopic Composition"]

    for i in 2:length(composition)
        elemNext = string(composition[i])
        massesNext = elementsAbundances[elemNext]["Relative Atomic Mass"]
        abundNext = elementsAbundances[elemNext]["Isotopic Composition"]

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
          masses = [eval(Meta.parse(charge_polarity))(masses[j], -1*elementsAbundances["e"]["Relative Atomic Mass"][1]) for j in 1:length(masses)]
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

isopat(x) = isotopicPattern(x)

isotopicPatternProtonated(x) = isotopicPattern(x, adduct="H", charge="+")
isopatprot(x) = isotopicPatternProtonated(x)

monoisotopicMass(x) = isotopicPattern(x;niceOutput=false)[1,1]
monoisomass(x) = monoisotopicMass(x)