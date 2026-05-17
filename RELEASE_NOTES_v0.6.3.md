# IsotopicCalc.jl v0.6.3 Release Notes

IsotopicCalc.jl v0.6.3 is a patch release focused on charge/adduct correctness and input validation.

## Fixed

### Isotopic pattern adduct handling

- `isotopic_pattern` now handles deprotonation/loss adducts such as `"H-"` consistently with `find_formula`.
- Radical and multiply charged adducts such as `"+2"` and `"2H+2"` now return charged m/z values instead of raising internal `BoundsError`s.
- Electron mass correction and charge-state division are now applied consistently for charged isotopic patterns.

Examples:

```julia
isotopic_pattern("H2O"; adduct="H-", print=false)    # [M-H]-
isotopic_pattern("C3H6O"; adduct="+2", print=false)  # [M]2+
isotopic_pattern("C3H6O"; adduct="2H+2", print=false) # [M+2H]2+
```

### Formula search with multiply charged ions

- `find_formula` now derives atom-count search bounds from the neutral mass implied by the observed m/z, charge state, and adduct mass.
- Valid formulas for multiply charged ions are no longer skipped when the neutral molecule is heavier than the measured m/z.

### Validation

- Invalid formulas such as `"2H"`, `"H0"`, `"[]"`, and `"(H)0"` now raise `ArgumentError` instead of lower-level parser or rounding errors.
- Invalid zero-charge adducts such as `"+0"` now raise `ArgumentError`.

## Tests

Added regression coverage for:

- Deprotonated isotopic patterns
- Radical and multiply charged isotopic patterns
- Multiply charged `find_formula` search bounds
- Malformed formula and adduct inputs

