# Release Notes for IsotopicCalc v0.6.2

## Bug Fixes

### Formula Finding - Low H/C Ratio Compounds

**Fixed**: `find_formula` now correctly identifies compounds with low hydrogen-to-carbon ratios, such as CO2, CO, and other small molecules.

#### Problem
Previously, the chemical heuristic rules in `find_formula` were too restrictive, requiring that formulas with carbon must have at least H ≥ 0.5×C. This eliminated many valid chemical compounds from the search results, including:
- Carbon dioxide (CO2)
- Carbon monoxide (CO)
- Other low H/C ratio compounds commonly found in mass spectrometry

For example, querying `find_formula(45.0; adduct="H+", tolerance_ppm=15000)` would miss CO2 (expected m/z 44.998, ~52 ppm error) even though it fell well within the tolerance range.

#### Solution
Modified Rule 1 in the chemical heuristic filtering to remove the lower bound requirement (H ≥ 0.5×C) while maintaining the upper bound (H ≤ 2×C+2+N). This allows formulas with H=0 while still pruning unrealistic high-hydrogen formulas.

**Technical Details:**
- **File**: `src/FindFormula.jl`
- **Change**: Simplified hydrogen ratio rule to only enforce upper bound
- **Impact**: More comprehensive formula search results, especially for small molecules and inorganic compounds
- **Performance**: Minimal impact; other heuristic rules still prune the search space effectively

#### Example
```julia
julia> find_formula(45.0; adduct="H+", tolerance_ppm=15000)
Matching formulas:
 CO2     [M+H]+  m/z: 44.997654  ppm: 52.15
 C2H4O   [M+H]+  m/z: 45.034040  ppm: -755.04
 CH2NO   [M+H]+  m/z: 45.020915  ppm: -464.57
 ...
```

**Before this fix**, CO2 and similar compounds would not appear in the results.

#### Testing
Added comprehensive regression tests to ensure:
- CO2 is detected at appropriate m/z values
- C2H4O and other valid formulas are still found
- All existing functionality remains intact
- Mass accuracy calculations are correct

## Changes

### Modified Files
- `src/FindFormula.jl` - Updated chemical heuristic Rule 1
- `test/runtests.jl` - Added "find_formula - Low H/C Ratio Compounds" test set

### Backward Compatibility
This is a **backward-compatible bug fix**. All existing code will continue to work as before, but users may see additional (previously missing) valid formulas in their search results.

## Migration Notes

No migration needed. If you were working around this limitation by using custom atom pools or post-processing results, you may now be able to simplify your code.

## Credits

This fix addresses user-reported issues where expected small molecules were missing from formula search results.

---

**Version**: 0.6.2
**Date**: 2025-12-13
**Type**: Bug Fix Release
