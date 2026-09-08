% Generated manually from R/capability-registry.R. Keep synchronized through
% tests/testthat/test-p0-capability-registry.R.
\name{longitudinal_capabilities}
\alias{longitudinal_capabilities}
\title{Supported longitudinal model capabilities}
\usage{
longitudinal_capabilities(component = c("routes", "margins", "copulas"))
}
\arguments{
\item{component}{Registry component to return: supported margin/copula
combinations (\code{"routes"}), marginal families (\code{"margins"}), or
copulas (\code{"copulas"}).}
}
\value{
A documentation-ready data frame. Every row contains the registry version so
exported tables retain their provenance.
}
\description{
Returns the versioned capability registry used by
\code{gamlss_longitudinal()} and the distribution-selection helpers. It
combines explicitly validated paper routes with dynamically discovered pure
continuous and pure discrete \pkg{gamlss.dist} families that satisfy the
package's generic family-interface contract.
}
\details{
Each fitted model uses one homogeneous marginal family across visits.
Continuous routes use a density-times-copula-density likelihood. Supported
count routes use exact discrete-discrete copula rectangle probabilities.
Bounded/binomial families are supported when the number of trials is supplied
explicitly with \code{margin_args = list(bd = ...)}. Ordinal/multinomial,
mixed continuous-discrete, logarithmic-series, and families that do not meet
the generic interface contract are rejected before optimization.

The public paper routes are BCPE margins with a Student-t copula and NBI
margins with a Clayton copula. The \code{limitations} and diagnostic columns
qualify the support available for other registered routes. Dynamically
discovered families use the Gaussian copula compatibility route; additional
copulas are listed only where that combination has explicit coverage.

The \code{hessian_path} column records an implemented analytical/numerical
computation path. It is not a guarantee of valid inference: the fitted Hessian
must pass the package's fit-level curvature validation. Evidence identifiers in
the route table name the tests that exercise each allowlisted route and the
generic diagnostic and Hessian machinery.
}
\examples{
longitudinal_capabilities("routes")
longitudinal_capabilities("margins")
longitudinal_capabilities("copulas")
}
