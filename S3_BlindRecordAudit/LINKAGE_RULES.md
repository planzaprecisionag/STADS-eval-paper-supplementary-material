# Blind record reconstruction and linkage rules

## Scope and authority

The final point-pair workbook defines the retained set and its presentation content. The delivered association GeoJSON records the candidate reconstruction. The rules below describe the supplied generate_blind_pair_association.py. The audit does not establish that every inferred association was independently verified in the field. Automated keyword-based agronomic labels in that script are not reviewer classifications and are not used as outcomes here.

## Reconstructing the candidate records

1. Read the corn and soybean prescribed-point layers and scouting exports. Group prescribed pairs by source layer, pair identifier, field, and farm. Pair numbers alone are not unique across the planned layers.
2. Collapse duplicate scouting exports using scouting-record identifier, field identifier, report identifier, comments, and the observation signature. The signature excludes picture URLs.
3. Extract valid nonzero numeric x/y coordinates from observation objects. The script does not itself read photograph EXIF metadata. These coordinates may have been populated upstream from photographs or other collection methods.
4. Limit candidate records to the same field identifier. Reject conflicting farm identifiers when both are available. Calculate the minimum haversine distance between the prescribed point and any coordinate in a candidate observation record.
5. Assign the nearest candidate within 35 m as a direct match. Record up to five nearby candidates within 75 m for inspection. This is an association tolerance, not an estimate of device accuracy or proof that an observation lies on the intended side of an anomaly boundary.

## Inference when a direct match is unavailable

Score candidates using the following additive evidence. A matching A/B label in comments receives +5; a conflicting label receives -2. Matching farm and field names receive +2 and +3, respectively. Relative to each directly matched partner record, matching scout, report, and submission timestamp receive +2, +4, and +2. An adjacent numeric record identifier receives +5; a difference of at most six receives +2 instead. The directly matched partner itself is excluded from the candidates for this step.

The highest-scoring candidate is accepted at a score of at least 7. This first inference pass has no spatial-distance requirement and can operate without a directly matched partner. A subsequent fallback uses a score of at least 8 and requires absent coordinates or a distance within 75 m. Candidates are evaluated in the source script's record ordering; tied scores do not establish a unique evidential match. The script does not impose one-to-one assignment across prescribed points. Reused observation records are flagged rather than automatically removed.

For inferred links, the delivered match_distance_m can describe a nearby spatial candidate rather than the selected inferred record. Consequently, the public direct_distance_m field is populated only for direct matches. Inference should not be interpreted as a measured positional error.

## Reconciliation to the final workbook

All 70 final presentations were matched to a unique prescribed point using pair identity and coordinates. They represent 55 unique pairs: 55 target-primary presentations and 15 additional reference-primary presentations. Only the 55 target-primary presentations enter the principal pair-level analyses. The other 15 presentations assess response consistency after reversing the observation columns.

All 14 inferred point links occur in the retained set. Primary observation descriptions were unchanged in the final workbook. One target role label was corrected from Int to In; the associated final presentation includes its paired observation comments. Author-confirmed intentional reuse is documented in data/blind_adjudication_log.csv. Other record reuse remains visible in data/blind_point_linkage_audit.csv. Pair and record identifiers in this package are pseudonyms; confidential crosswalks and coordinates are not distributed.

The public analysis script reproduces outcomes and audit summaries from the delivered de-identified associations. Re-running spatial linkage itself requires the confidential raw coordinates and scouting records; it cannot be done from the public audit alone.
