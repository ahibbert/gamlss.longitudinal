# Frozen review bundle

## Source identity

- Canonical manuscript repository: `ahibbert/working-paper-gamlss-long`
- Canonical manuscript commit: `68c3bad26626ce7c267bd330364cfda8df7a6b76`
- Package checkout commit: `7343d06c986b81e2dbbbf5c2bc2732c3fe539c41`
- Package checkout state: dirty; the review intentionally includes the current
  uncommitted replication and package work described by the staged assets.
- Review date: 2026-08-27, Australia/Sydney.

## Build process

1. A fresh depth-one clone of the canonical manuscript repository was made.
2. `Rscript paper/replicate.R --profile paper --workers 1 --no-restore --no-install`
   refreshed the local public paper outputs.
3. `paper/publish-assets.R` was run first with `--dry-run` and then with
   `--apply` against the disposable clone. Only manifest-allowlisted public
   `charts/` and `tables/` assets were staged.
4. Three previously documented mechanical compilation fixes were applied only
   to the disposable clone: three `\textt` misspellings were changed to
   `\texttt`, a paragraph break was removed from a revision-note argument, and
   the second duplicate `Sklar1973` bibliography key was renamed.
5. `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` produced the
   frozen 47-page PDF.

No substantive manuscript wording was changed, and the live Overleaf repository
was not modified.

## Frozen artifact hashes

| Artifact | SHA-256 | Size |
|---|---|---:|
| Compiled `main.pdf` | `3c4e4cfe08e65cd79dcd8e95581d933267f62716de4165b5bc00cfc45f6609b9` | 2,848,504 bytes |
| Compiled-clone `main.tex` | `b34dfe7e8ed13b6e12a1868d5b9ff330d4e6a317257dafb3d04fcfcfb27648ee` | - |
| Compiled-clone `sample.bib` | `4f6e5d0385e9ff14f46d87ce178e380aca616bc1a49c9d898521db09d404a02d` | - |
| Public paper manifest | `22343b7532bc0e108036887c6c0eb131f6337275565da0174132c79904c23f4b` | - |

The frozen build is located at
`C:/Users/Aydin/AppData/Local/Temp/gamlss-jss-review-round-20260827/paper-repo/main.pdf`.
All 47 pages were rendered to PNG for visual review.

## Known mechanical debt retained for reviewers

- Eight unresolved citation occurrences involving six keys: `czado2015`,
  `Beareseo2015`, `lambert_copula`, `topmodels`, and `marra2025?`.
- Seventeen overfull hboxes, 12 underfull hboxes, and 12 underfull vboxes.
- Red revision notes and author questions remain visible because this is a
  developmental review rather than a submission simulation.

## Review independence

Each specialist reviewer received the same PDF/TeX snapshot and a role-specific
brief. Reviewers were instructed not to read existing internal manuscript audits
or other reviewers' reports. The final synthesis may use current official JSS
author and style guidance and may inspect the repository to verify findings.
