# Online Resource S3

Blind-record reconstruction, linkage audit, secondary statistical methods, and sensitivity analyses for the STADS field-evaluation case study.

Start with **S3_Narrative.docx**. The package contains de-identified analysis inputs, machine-readable results, documented linkage rules, and a standalone reproduction script. Current RTK analysis records are included to reproduce the stream-specific secondary summaries; the full RTK reconciliation and correction audit belongs in **Online Resource S8**.

## Files and roles

- `S3_Narrative.docx`: publication-oriented methods and results narrative.
- `LINKAGE_RULES.md`: reconstruction rules, inference thresholds, reconciliation decisions, and limits of public reproducibility.
- `DATA_DICTIONARY.md`: units, identifiers, missing values, endpoints, and table conventions.
- `data/blind_pair_inventory.csv`: all 99 prescribed pairs and retained-versus-excluded attributes.
- `data/blind_point_linkage_audit.csv`: all 198 prescribed points, link methods, observation reuse, and direct-match distances.
- `data/blind_inferred_link_audit.csv`: the 14 inferred links, extracted for inspection.
- `data/blind_adjudication_log.csv`: role correction, deliberate reuse, and presentation-key reconciliation.
- `data/blind_reviewer_ratings.csv`: all 350 classifications, including the 75 ratings from 15 reversed presentations.
- `data/reviewer_subsets.csv`: the competing-interest and historical masking subsets, which differ.
- `data/rtk_analysis_records.csv`: the current 57 retained RTK pairs.
- `data/source_provenance.csv`: source-version fingerprints without confidential file identities.
- `results/`: regenerated consensus, sensitivity, agreement, confidence, retention, reversal, and validation tables.
- `results/rtk_endpoint_summary.csv`: authoritative supplied RTK summary; all 11 rows are checked in `rtk_summary_verification.csv`.
- `results/supplied_blind_bootstrap_intervals.csv`: blind-only rows from the supplied bootstrap output. Old RTK rows have been removed. Original Wilson comparison columns are retained and explicitly named; they are not the primary intervals.
- `results/supplied_reversal/`: sanitized historical reversed-presentation outputs. These preserve their original analysis definitions and are distinguished from regenerated summaries.
- `code/reproduce_s3.py` and `code/analysis_functions.py`: public analysis and extracted v5 reliability functions.
- `code/requirements.txt`: versions used to regenerate this package.
- `MANIFEST_SHA256.csv`: integrity manifest for the files in this resource, excluding the manifest itself.

## Reproduce

With Python 3.12 and the dependencies in code/requirements.txt, run from this directory:

```bash
python -m pip install -r code/requirements.txt
python code/reproduce_s3.py
```

The script overwrites the regenerated CSVs in results/ and writes validation.json. It does not overwrite the supplied historical files. Bootstrap seeds are explicit per row. Some secondary intervals were generated for this package; they should not be represented as previously archived analyses. Complete cases mean all five reviewers supplied an evaluable presence classification.

## Interpretation

The two scouting streams use different samples and procedures. Their percentages are complementary evidence, not a randomized or matched comparison of scouting performance. Endpoint-specific exclusions preserve usable evidence for other endpoints. Reversed presentations and deliberately reused records are not additional independent field observations. Field-cluster resampling retains dependence within fields but does not remove selection bias, dependence between fields, or uncertainty in inferred links.

## Other online resources

S1 provides blind observation content; S2 provides the RTK collection form and observations; S4 and S5 provide the blank reviewer workbook and rubric; S6 provides completed reviewer workbooks and the broader source processing scripts; S7 provides field examples; S8 provides RTK reconciliation. S3 avoids duplicating private raw exports or obsolete analysis versions.
