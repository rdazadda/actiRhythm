## Submission

This is a new submission. actiRhythm is a standalone package for circadian
rest-activity rhythm analysis of activity time series (it is the analysis engine
behind a larger actigraphy dashboard, released here as an independent package).

## Test environments

- Windows 11, R 4.5.3 (local)
- GitHub Actions: ubuntu-latest, windows-latest, and macos-latest, on R release
  and R-devel
- win-builder, R-devel

## R CMD check results

Local `R CMD check --as-cran` (Windows 11, R 4.5.3): 0 errors | 0 warnings | 1 note.

The one note is the expected new-submission note:

```
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Raymond Dacosta Azadda <rdazadda@alaska.edu>'
  New submission
```

The `URL` and `BugReports` fields point to the package's public GitHub
repository, which will be live at the time of submission. Any words flagged by
the spell check are correctly spelled chronobiology terms or author surnames from
the method references, all listed in `inst/WORDLIST`.

## Notes

- No `\dontrun{}` is used. Examples either run quickly or, where they read the
  bundled sample data or perform a full multi-method analysis, are wrapped in
  `\donttest{}`.
- The package implements published circadian methods, cited by author and year
  in the Description and in full in each function's `\references` section.
- Functions that draw bootstrap samples take an optional `seed` argument passed
  to `set.seed()`; they do not otherwise alter the global random-number state.
