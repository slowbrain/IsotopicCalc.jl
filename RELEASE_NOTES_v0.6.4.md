# IsotopicCalc.jl v0.6.4 Release Notes

IsotopicCalc.jl v0.6.4 is a patch release focused on stricter validation for isotope formula syntax.

## Fixed

### Square-bracket isotope validation

- Malformed square-bracket isotope syntax now raises `ArgumentError` consistently.
- Stray closing brackets such as `"C]H4"` are no longer silently ignored.
- Nested or malformed isotope tokens such as `"[[13C]]H4"`, `"[C13]H4"`, `"[13]H4"`, and `"[13c]H4"` are rejected before mass lookup.
- The non-isotope formula parser now rejects unexpected leftover syntax rather than skipping unrecognized characters.

## Tests

Added focused regression coverage for:

- Stray closing square brackets
- Doubled closing square brackets
- Nested square brackets
- Malformed isotope token contents and casing

