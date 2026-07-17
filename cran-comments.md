## Submission

This is a new submission. actiRhythm quantifies the circadian rest-activity
rhythm from activity counts and raw accelerometer recordings.

## Test environments

- win-builder, R-devel (R Under development, 2026-07-16 r90264 ucrt),
  Windows Server 2022
- Windows 11, R 4.5.3 (local), `R CMD check --as-cran`
- GitHub Actions: ubuntu-latest, windows-latest, and macos-latest, on R release
  and R-devel

## R CMD check results

win-builder R-devel (R Under development, 2026-07-16 r90264 ucrt): 0 errors |
0 warnings | 1 note. The PDF and HTML manuals, the vignettes, the examples, and
the test suite all build and pass.

The one note is the expected new-submission note:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Raymond Dacosta Azadda <rdazadda@alaska.edu>'

New submission

Possibly misspelled words in DESCRIPTION:
  Bushell (24:32)
  Cornelissen (21:50)
  ENMO (33:61)
  Hees (32:52, 34:65)
  Kantelhardt (27:5)
  Lomb (23:51)
  Marler (22:43)
  Peng (26:29)
  Scargle (23:56)
  Sokolove (24:19)
  Someren (19:19)
  Winnebeck (30:19)
  al (16:73, 19:30, ...)
  cosinor (20:5, 21:5)
  et (16:70, 19:27, ...)
  interdaily (15:66)
  intradaily (16:24)
```

The flagged words are all spelled correctly: author surnames from the method
references (Bushell, Cornelissen, Hees, Kantelhardt, Lomb, Marler, Peng,
Scargle, Sokolove, Someren, Winnebeck), standard chronobiology terms (cosinor,
interdaily, intradaily, ENMO), and "et al." from the citations. They are listed
in `inst/WORDLIST`.

The `URL` and `BugReports` fields point to the package's public GitHub
repository, which is live.

## Notes

- No `\dontrun{}` is used. Examples either run quickly or, where they read the
  bundled sample data or perform a full multi-method analysis, are wrapped in
  `\donttest{}`.
- The package implements published circadian methods, cited by author and year
  in the Description and in full in each function's `\references` section.
- Functions that draw bootstrap samples take an optional `seed` argument passed
  to `set.seed()`; they do not otherwise alter the global random-number state.
