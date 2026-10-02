# Data dictionary

## Conventions

CSV files are UTF-8 with one header row. True/False values are Boolean indicators. Blank numeric or categorical cells mean missing or not evaluable, not a negative result. Percent fields use 0–100 units; alpha, kappa and AC1 are unitless. Distances are meters. Dates describe imagery scenes, not necessarily visits. File-specific columns below define the main inputs and regenerated outputs.

## Blind identifiers and inventory

`pair_id` (BL-): pseudonymous prescribed pair, unique across all source layers. `field_id` (BF-) and `farm_id` (FARM-) preserve grouping without revealing names. `observation_id` (OBS-) identifies a linked source scouting record and makes reuse visible; blank means no linked record. These pseudonyms are local to S3 and are not a public crosswalk to S1 or S6. RTK identifiers retain the supplied S8 public identifiers.

`crop_source`: corn/soybean prescribed source layer. This differs from `crop`, the crop recorded in the final reviewer workbook, which can include other crops. `window_id`: source sampling-window code. `scene_date`: source imagery date. `stads_direction`: mapped direction at the target point, normalized from plus/minus or positive/negative labels. `planned_separation_m`: haversine distance between prescribed target and reference. `retained`: pair appears in final authoritative workbook. `linked_points`: number of associated points. `inferred_points`: number associated by contextual inference. `record_reuse`: at least one source record used at multiple points. `retention_reason`: final-workbook inclusion status; it is not a newly adjudicated reason for missing photographs.

## Point audit and ratings

`point_role`: target (inside) or reference (outside), including the confirmed Int correction. `point_letter`: original A/B label, which is not always unique within a pair. `match_method`: direct, inferred_adjacent_or_point_label, or unmatched. `direct_distance_m`: source minimum observation-coordinate distance for direct links only. `observation_coordinate_count`: coordinate-bearing observations in selected record. `record_reuse_count`: delivered reuse count, not an independent-observation count. `role_typo_corrected`: confirmed Int-to-In correction.

`reviewer_id`: R1–R5, aligned to the supplied reversal outputs. `presentation`: role of the primary column shown to the reviewer. Target-primary ratings are used for principal analyses. `presence`: Confirmed, Possibly Confirmed, Not Confirmed, Indeterminate, or Excluded. Only the first three are evaluable. `direction`: Positive, Negative, Same, Indeterminate, or Excluded. Positive/Negative are evaluable for binary mapped-direction agreement. Same can contribute to the direction mode but cannot agree with a binary STADS direction. `presence_confidence` and `direction_confidence`: Low, Medium, High, or unavailable/indeterminate.

`identified_competing_interest`: the two author-identified reviewers excluded in the three-reviewer sensitivity analysis. `historical_masking_subset`: designation preserved from supplied reversal outputs; it does not independently establish a reviewer's exposure to maps.

## Consensus and endpoints

`n_presence_valid` and `n_direction_valid`: counts of evaluable reviewer responses. `presence_evaluable`: at least two valid presence responses. `Strict`: at least ceil(0.75*n) Confirmed. `Moderate`: at least ceil(0.75*n) Confirmed or Possibly Confirmed. `Liberal`: at least ceil(0.50*n) Confirmed or Possibly Confirmed. `confidence_weighted_support`: weighted mean of support scores 1, 0.5, 0 with confidence weights 1, 2/3, 1/3. `Confidence-weighted`: support at least 0.5 with at least two valid presence ratings. Valid confidence is required for a rating to contribute to the weighted mean.

`presence_mode`/`consensus_direction`: unique most frequent evaluable category; Tie denotes multiple maxima. `direction_evaluable`: unique binary consensus mode. `direction_match`: unique consensus equals STADS direction. `composite`: Full correspondence when Moderate presence and direction agree; Presence-only correspondence when Moderate presence holds but direction does not agree or cannot be determined; No correspondence when presence is evaluable and Moderate fails. Two presence-nonevaluable pairs have blank composite labels.

`subset`, `endpoint`, and `group` define each summary row. `n` is the endpoint denominator, `successes` the numerator, `percent` their ratio, `n_fields` the number of eligible field clusters. `ci_lower_percent`/`ci_upper_percent` are 95% percentile field-cluster bootstrap bounds. `seed` is the NumPy default_rng seed. No interval is estimated with fewer than two fields. Crop summaries are descriptive, including very small strata.

## RTK

`public_pair_id` and `public_field_cluster_id` retain supplied public identifiers. `inside_comparison`/`outside_comparison` contain Better, Worse, or Same comparisons. The *_evaluable flags identify endpoint-specific eligibility; `anomaly_confirmed`, `direction_agreement`, and `composite_classification` are the current supplied outcomes. `observed_direction` is the interpreted field direction and `stads_direction` is the mapped direction. RTK outcome rules and exclusions are described in the narrative and in S8.

## Reliability and confidence

Reliability rows specify the category set. `items_two_or_more` is the number available for coincidence-based coefficients. Fleiss' kappa uses the modal number of nonmissing ratings; `fleiss_items` reports that restricted denominator. Pairwise Cohen's kappa uses each reviewer pair's common items. Raw agreement is a 0–1 fraction. The interval alpha uses squared distances among the presence support scores. Confidence summaries report support, agreement with the panel mode including the same reviewer, or binary direction agreement. These are descriptive associations, not external calibration.

## Reversal and retained historical tables

Regenerated reversal tables compare target-primary with reference-primary presentations. Positive and Negative are swapped before direction comparison; Same, Indeterminate and Excluded remain unchanged. All-category agreement therefore includes agreement on indeterminate responses. The historical tests retain their original labels and confidence intervals, including non-cluster-adjusted intervals. They are secondary descriptive checks and should not replace the field-cluster intervals for primary endpoints. `my_pair_id` was replaced with `pair_id`; free-text rationales, source row identifiers, and coordinates were removed. No rationale flags or detection warnings were present in the supplied files; their absence is not proof of blinding.
