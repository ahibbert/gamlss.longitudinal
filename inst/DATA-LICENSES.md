# Data licenses and attribution

The `gamlss.longitudinal` package is GPL-3. The three analysis-ready datasets
below are derived from third-party public datasets and retain the stated
source licenses and attribution. The preparation script is
`data-raw/build-public-candidate-data.R` in the source repository.

## `patents_panel`

- Upstream object: `pglm::PatentsRDUS`, pglm 0.2-4.
- Upstream package license: GPL >= 2.
- Source-package SHA-256:
  `caf085e7f5d693efdb06a96057c6e1d9e132a9525b0bd6a595804e2b853659e3`.
- Study: Hall, B. H., Griliches, Z., and Hausman, J. A. (1986), "Patents
  and R&D: Is There a Lag?", *International Economic Review*, 27, 265--283.
- Documentation: <https://stat.ethz.ch/CRAN/web/packages/pglm/refman/pglm.html>.
- Modifications: reshaped/selected into an analysis-ready long panel, original
  CUSIP identifiers replaced by sequential integers, variable names clarified,
  and the outcome-derived `sumpat` field omitted.

## `pbc_prothrombin`

- Upstream object: `survival::pbcseq`, survival 3.8-11.
- Upstream package license: LGPL >= 2.
- Source-package SHA-256:
  `4a87aea323d477e142c36601509f3771f150b441df6db75537ea1777e6546888`.
- Copyright: 2000 Mayo Foundation for Medical Education and Research.
- Documentation: <https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pbcseq.html>.
- Study reference: Murtaugh, Dickson, Van Dam, Malinchoc, Grambsch,
  Langworthy, and Gips (1994), "Primary biliary cirrhosis: prediction of
  short-term survival based on repeated patient visits", *Hepatology* 20,
  126--134.
- Modifications: selected modelling variables, patient identifiers normalized
  to sequential integers, baseline variables repeated by participant,
  visit order added, and elapsed years calculated from follow-up days.

The GNU Lesser General Public License is available at
<https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.

## `vietnam_steps`

- Authors: Tram T. N. Truong, Van-Anh N. Huynh, and Kien G. To.
- Article: "Comparison between self-reported and pedometer-measured physical
  activity in Vietnamese adolescents: A reliability and agreement study"
  (2025), *PLOS Global Public Health*.
- Article DOI: <https://doi.org/10.1371/journal.pgph.0004725>.
- Dataset DOI: <https://doi.org/10.1371/journal.pgph.0004725.s002>.
- License: Creative Commons Attribution 4.0 International (CC BY 4.0),
  <https://creativecommons.org/licenses/by/4.0/>.
- Source-file MD5: `722cb5aa627ff9273712c485e130800b`.
- Modifications: converted to long format; original study identifiers removed;
  BMI and PAQ-C derived from the documented source fields; two step counts
  outside the article's valid range and three visibly fractional counts changed
  to `NA`; source-cell status retained in `step_status`; unused raw fields
  omitted.
