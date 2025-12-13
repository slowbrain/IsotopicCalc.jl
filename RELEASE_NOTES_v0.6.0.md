# IsotopicCalc.jl v0.6.0 Release Notes

We're excited to announce **IsotopicCalc.jl v0.6.0**, featuring a major enhancement to the adduct system with flexible parsing and multiple charge support!

## 🎯 Highlights

### 1. Flexible Adduct System - Use Any Element! 🧪

The adduct system has been completely redesigned to support **any element from the periodic table**, not just the hardcoded H, Na, and K.

**Before (v0.5.0):**
```julia
# Only supported: H+, Na+, K+, H-, +, ""
find_formula(59.0491; adduct="H+")
find_formula(81.0284; adduct="Na+")
```

**Now (v0.6.0):**
```julia
# Use ANY element!
find_formula(99.0; adduct="Ca+")   # Calcium
find_formula(65.0; adduct="Li+")   # Lithium
find_formula(82.0; adduct="Mg+")   # Magnesium
find_formula(94.0; adduct="Cl-")   # Chloride loss
find_formula(100.0; adduct="Br-")  # Bromide loss
```

### 2. Multiple Charge Support ⚡⚡⚡

Critical for ESI-MS and peptide analysis! Now you can specify multiply charged ions.

**New Notation:** `[n]Element+/-[m]`
- `n` = number of atoms (optional, defaults to 1)
- `m` = charge state (optional, defaults to 1)

**Examples:**
```julia
# Doubly protonated (common in ESI-MS)
find_formula(30.0; adduct="2H+2")
# Output: C3H6O  [M+2H]2+  m/z: 30.028037

# Triply protonated (peptides)
find_formula(20.0; adduct="3H+3")
# Output: C2H3NO  [M+3H]3+  m/z: 20.014431

# Radical dication
find_formula(29.0; adduct="+2")
# Output: C2H2O2  [M]2+  m/z: 29.002191

# Radical trication
find_formula(19.5; adduct="+3")
```

### 3. Important Distinction: Atoms vs. Charge

Understanding the difference between atom count and charge state:

```julia
# "2H+" = 2 protons added, charge state +1
find_formula(60.0; adduct="2H+")
# m/z = (M + 2×1.008) / 1 = M + 2.016

# "2H+2" = 2 protons added, charge state +2
find_formula(30.0; adduct="2H+2")
# m/z = (M + 2×1.008) / 2 = (M + 2.016) / 2
```

### 4. Type Flexibility

No more type errors with integer inputs!

```julia
# Both work now
find_formula(33; adduct="H+")      # Integer ✓
find_formula(33.0; adduct="H+")    # Float ✓
```

## 📋 Complete Adduct Notation Reference

### Format
```
[n]Element+/-[m]
```

### Single Charge Examples
| Notation | Description | Display |
|----------|-------------|---------|
| `"H+"` | Protonated | `[M+H]+` |
| `"Na+"` | Sodium adduct | `[M+Na]+` |
| `"K+"` | Potassium adduct | `[M+K]+` |
| `"Ca+"` | Calcium adduct | `[M+Ca]+` |
| `"2H+"` | Two protons, charge +1 | `[M+2H]+` |
| `"H-"` | Deprotonated | `[M-H]-` |
| `"Cl-"` | Chloride loss | `[M-Cl]-` |
| `"+"` | Radical cation | `[M]+` |
| `""` | Neutral | (no adduct) |

### Multiple Charge Examples
| Notation | Description | Display |
|----------|-------------|---------|
| `"2H+2"` | Doubly protonated | `[M+2H]2+` |
| `"3H+3"` | Triply protonated | `[M+3H]3+` |
| `"4H+4"` | Quadruply protonated | `[M+4H]4+` |
| `"+2"` | Radical dication | `[M]2+` |
| `"+3"` | Radical trication | `[M]3+` |
| `"-2"` | Radical dianion | `[M]2-` |
| `"2H-2"` | Double deprotonation | `[M-2H]2-` |

## 🔬 Real-World Applications

### Peptide Analysis
```julia
# Typical peptide with MW ~600 Da
# Doubly charged: m/z ≈ 301
find_formula(301.0; adduct="2H+2", atom_pool=Dict(
    "C"=>50, "H"=>100, "N"=>20, "O"=>20, "S"=>2
))

# Triply charged: m/z ≈ 201
find_formula(201.0; adduct="3H+3", atom_pool=Dict(
    "C"=>50, "H"=>100, "N"=>20, "O"=>20, "S"=>2
))
```

### Metal Adducts
```julia
# Lithium adducts (common in MALDI)
find_formula(65.0; adduct="Li+")

# Calcium adducts (biological samples)
find_formula(99.0; adduct="Ca+")

# Magnesium adducts
find_formula(82.0; adduct="Mg+")
```

### Halogen Loss
```julia
# Chlorine loss (organochlorines)
find_formula(94.0; adduct="Cl-")

# Bromine loss (organobromines)
find_formula(100.0; adduct="Br-")
```

## 🔧 Technical Details

### Implementation
- **Dynamic parsing**: Replaced hardcoded `ADDUCTS` dictionary with `parse_adduct()` function
- **Element validation**: Automatically validates against all elements in the periodic table
- **Regex pattern matching**: `^(\d*)([A-Z][a-z]?)([+-])(\d*)$` for flexible parsing
- **Backward compatible**: Old `get_adduct_info()` function still works via wrapper

### m/z Calculation
The m/z calculation correctly handles multiple charges:

```julia
M = neutral_mass + adduct_mass - charge × 0.0005485  # electron mass correction
m/z = charge == 0 ? M : M / abs(charge)
```

### Display Formatting
- Single charge: `[M+H]+`, `[M+Na]+`
- Multiple charges: `[M+2H]2+`, `[M+3H]3+`, `[M]2+`
- Proper superscript notation in output

## 📦 Upgrading

### From v0.5.0
No breaking changes! All existing code continues to work:

```julia
# Still works exactly as before
find_formula(59.0491; adduct="H+")
find_formula(81.0284; adduct="Na+")
find_formula(58.0419; adduct="+")
```

### Installation
```julia
using Pkg
Pkg.update("IsotopicCalc")
```

## 🐛 Bug Fixes
- Fixed JSON.Object type compatibility issue with SubString keys
- Improved error messages for invalid adduct formats

## 📚 Documentation
- Comprehensive docstrings with examples
- Updated README with new notation
- Added CHANGELOG entry with migration guide

## 🙏 Acknowledgments

This release brings IsotopicCalc.jl's adduct system to a new level of flexibility, making it suitable for a much wider range of mass spectrometry applications including ESI-MS, MALDI-TOF, and peptide/protein analysis.

---

**Full Changelog**: [v0.5.0...v0.6.0](https://github.com/slowbrain/IsotopicCalc.jl/compare/v0.5.0...v0.6.0)

**Questions or Issues?** Please open an issue on [GitHub](https://github.com/slowbrain/IsotopicCalc.jl/issues)
