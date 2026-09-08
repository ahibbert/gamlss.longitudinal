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
combines explicitly validated paper routes with a generic compatibility
contract for pure continuous and pure discrete \pkg{gamlss.dist} families.
}
\details{
Each fitted model uses one homogeneous marginal family across visits.
Continuous routes use a density-times-copula-density likelihood. Discrete
routes use exact discrete-discrete copula rectangle probabilities. Binomial
families require explicit trial counts through
\code{margin_args = list(bd = ...)}. Ordinal/multinomial,
logarithmic-series, and mixed continuous-discrete families are rejected before
optimization because they require different response or likelihood machinery.

The public paper routes are BCPE margins with a Student-t copula and NBI
margins with a Clayton copula. The \code{limitations} and diagnostic columns
qualify the support available for other registered routes.

The \code{hessian_path} column records an implemented analytical/numerical
computation path. It is not a guarantee of valid inference: the fitted Hessian
must pass the package's fit-level curvature validation. Evidence identifiers in
the route table distinguish explicitly validated routes from generic family
interface and bounded-fit probes.
}
\examples{
longitudinal_capabilities("routes")
longitudinal_capabilities("margins")
longitudinal_capabilities("copulas")
}
