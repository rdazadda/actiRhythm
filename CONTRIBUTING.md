# Contributing to actiRhythm

Thanks for your interest in actiRhythm. A few notes for contributors.

## Questions and bugs

Open an issue at <https://github.com/rdazadda/actiRhythm/issues>. For a bug, a small
reproducible example, a short count vector with timestamps or one of the bundled
`example_agd()` files, makes it much faster to diagnose.

## Pull requests

- Branch from `main` and keep each pull request focused on one change.
- Run `devtools::check()` and confirm the test suite passes before opening it.
- Match the surrounding code style, and add a test for any new behaviour.
- A new circadian method should cite the paper that defines it in the function's
  `@references`.

## Conduct

Please be respectful in all project spaces.
