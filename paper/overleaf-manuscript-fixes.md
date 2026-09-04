# Overleaf manuscript fixes found during replication QA

The asset publisher does not edit manuscript source. Apply these source-only fixes
manually in the canonical Overleaf repository before the final reviewer compile:

1. In the paragraph discussing sandwich inference, replace the three misspelled
   commands `\textt{summary()}`, `\textt{confint()}`, and
   `\textt{wald\_test()}` with `\texttt{summary()}`, `\texttt{confint()}`, and
   `\texttt{wald\_test()}`.
2. In the `\jssrevisionnote{...}` immediately below
   `\subsubsection{Mis-specification of correlation structure and benchmark to GEE}`,
   remove the blank paragraph between the sentences ending “covariate dependent
   correlation.” and beginning “I guess...”. A paragraph break is not valid inside
   that command argument and currently causes a runaway-argument compile failure.
3. Remove or rename one of the two `@article{Sklar1973, ...}` entries in
   `sample.bib`. The duplicate BibTeX key currently makes `bibtex main` exit with
   an error.
4. Add or correct bibliography entries for the still-undefined citation keys
   `czado2015`, `Beareseo2015`, `lambert_copula`, `topmodels`, and
   `marra2025?`. The compile succeeds after fixes 1--3, but reports eight
   unresolved citation occurrences involving these six keys.

These corrections were tested only in a disposable compile clone. They are not
made by `paper/publish-assets.R` and have not been applied to the live paper.
