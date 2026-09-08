# Missingness publication inputs

This directory is the portable, reviewer-facing output of the registered
missingness study. It contains 240 simulated tasks and 480 model attempts.
Subject-level monotone dropout is the headline analysis; time-dependent
intermittent MAR is labelled as a sensitivity analysis.

The CSVs report attempted, converged, retained, and failed fits; metric-specific
denominators; failure-inclusive sensitivity estimates; Monte Carlo uncertainty;
and the exact expected-warning audit. `missingness_configuration.csv` records
the fixed design, seeds, optimizer limits, and failure penalty.

Raw checkpoint payloads, machine paths, process identifiers, quarantine logs,
and signing machinery are deliberately not publication inputs. They remain
local producer infrastructure and are not needed to reproduce the manuscript
figures from these CSVs.
