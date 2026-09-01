# Legacy missingness bundle

These tracked files are retained only as historical review inputs. They use the
obsolete `mar` / `time_mar` study, omit genuine monotone dropout, and do not
contain the required attempt/failure, failure-inclusive sensitivity, or Monte
Carlo uncertainty artifacts.

The paper-profile loader deliberately rejects this directory. Replace the
bundle only with outputs from the current full missingness producer, including:

- `missingness_design_registry.csv`
- `missingness_estimand_registry.csv`
- `missingness_sensitivity_registry.csv`
- `missingness_checkpoint_status.csv`
- `missingness_checkpoint_payloads.rds`
- `missingness_checkpoint_content_manifest.csv`
- `fit_run_log.csv`
- `attempt_failure_summary.csv`
- `failure_reason_summary.csv`
- `missingness_by_rep.csv`
- `missingness_pattern_by_subject_visit.csv`
- `fixed_effects_by_rep.csv`
- `smooth_estimates_by_rep.csv`
- `missingness_headline_summary.csv`
- `fixed_term_summary_by_missingness.csv`
- `smooth_irmse_summary.csv`
- `smooth_selected_plot_data.csv`

Until that production run and promotion are complete, these CSVs must not be
used to support current manuscript claims or figures.
