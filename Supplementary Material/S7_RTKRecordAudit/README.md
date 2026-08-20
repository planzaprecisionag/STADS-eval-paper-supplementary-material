# Online Resource S7: RTK record reconciliation and analysis audit

This directory provides a de-identified, reproducible audit trail for the
RTK-GNSS-guided scouting records used in the STADS field-evaluation analysis.
It reconciles the original planned pairs, one replacement, completed Survey123
pairs, post-scouting exclusions, retained analysis pairs, endpoint-specific
denominators, and field-cluster bootstrap summaries.

## Record flow

- Original planned pairs: 85
- Replacement pairs added: 1
- Truly unscouted original pairs: 27
- Completed pairs: 59
- Completed non-target pairs excluded: 2
- Retained analysis pairs: 57
- Anomaly-confirmation/composite evaluable: 53
- Direction evaluable: 43

The files reproduce 44 confirmed anomalies,
37 direction matches, and composite counts of
37 full, 7
presence-only, and 9 no correspondence.

## Public files

- `rtk_record_reconciliation.csv`: one row for each of the 85 original planned
  pairs and the replacement pair.
- `rtk_analysis_records.csv`: one row for each of the 57 retained analysis
  pairs, with anonymized field-cluster IDs and categorical endpoint inputs and
  outcomes.
- `rtk_endpoint_summary.csv`: endpoint counts, percentages, cluster counts,
  field-cluster bootstrap confidence intervals, replicate counts, and seeds.
- `rtk_reconciliation_audit.json`: machine-readable source fingerprints,
  reconciled counts, and validation results.
- `rtk_reconciliation_audit.txt`: human-readable audit summary.

## Endpoint logic

Each retained pair contains one observation inside and one outside the mapped
anomaly. A pair is non-evaluable for anomaly confirmation when one member is
rated `Same` and the other is rated `Better` or `Worse`; this occurred for four
pairs. Otherwise, a pair is confirmed when at least one comparison indicates
`Better` or `Worse`, and it is not confirmed when the available comparisons
indicate only `Same`.

Observed direction is positive when the inside observation is `Better` or the
outside observation is `Worse`, and negative when the inside observation is
`Worse` or the outside observation is `Better`. Conflicting directional cues
are indeterminate; this occurred for one confirmed pair. Full correspondence
requires confirmed anomaly presence and direction agreement. Presence-only
correspondence indicates confirmed presence with mismatched or indeterminate
direction. No correspondence indicates that anomaly presence was not confirmed.

## Bootstrap reproduction

The public analysis file contains anonymized field-cluster identifiers. For
each endpoint/group, field clusters are sampled with replacement, all records
from each sampled cluster are retained, the percentage is recalculated for
10,000 replicates, and the 2.5th and 97.5th percentiles form the confidence
interval. The seeds are reported in `rtk_endpoint_summary.csv`.

## Privacy

The public files omit geometry, coordinates, anomaly polygons, original pair,
farm, region, and field identifiers, Survey123 GUIDs, creator/editor fields,
photograph and local file paths, free-text observations, correction notes, and
numeric anomaly scores. Public pair and field-cluster IDs have no released
crosswalk to internal identifiers.


