#%%
"""
STADS reviewer assessment analysis
=================================

This script extends the current reviewer-assessment aggregation workflow to run
six publication-oriented analyses:

1. Inter-rater reliability
   - Krippendorff's alpha for anomaly confirmation and direction labels
   - Fleiss' kappa where the rater matrix is sufficiently complete
   - Gwet's AC1 for nominal agreement
   - Pairwise Cohen's kappa among reviewer spreadsheets

2. Consensus-based STADS support
   - Strict, moderate, liberal, and confidence-weighted confirmation rules
   - Wilson confidence intervals for confirmation rates

3. Direction agreement
   - Reviewer-level direction agreement with STADS direction
   - Consensus-level direction agreement with STADS direction
   - Wilson confidence intervals overall and by crop / STADS direction

4. Mixed-effects / clustered logistic model
   - Binary confirmation outcome at the reviewer-rating level
   - Reviewer and field random effects if statsmodels BinomialBayesMixedGLM is available
   - GLM with reviewer fixed effects and clustered SE fallback

5. Blind vs RTK-guided comparison
   - Treats the four reviewer spreadsheets as the blind evidence stream
   - Treats the Survey123 GeoPackage as an RTK-guided, single-observer evidence stream
   - Uses sample_pai as the point-pair identifier and a_or_b as the within-pair A/B label
   - Infers RTK anomaly-point evidence from compare_to within each sample_pai pair

6. Reviewer confidence analysis
   - Confirmation and consensus agreement by confidence class
   - Direction agreement by direction-confidence class

Notes
-----
- The script preserves identical rows across spreadsheets as independent reviewer
  classifications.
- Point filtering follows the existing logic: keep only rows whose
  my_pair_id || my_pair_a_or_b composite key is labeled as in-anomaly in the
  location/evaluation CSV.
- Configure paths and column names in the CONFIG section below before running.
- RTK-guided observations are not included in inter-rater reliability because they
  represent one field observer rather than independent retrospective reviewers.

Dependencies
------------
Required: pandas, numpy, matplotlib, openpyxl
Optional: statsmodels for model fitting
"""

# from __future__ import annotations

import itertools
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# CONFIG
# =============================================================================

# Location/evaluation CSV used to determine whether each reviewed point is inside
# a STADS anomaly and to attach STADS direction.
FILEPATH_FOR_LOCATIONS = Path(
    r"<YOUR PATH HERE>/S4.stads_scouting_validation_analysis_all_point_data_with_evaluations.csv"
)

# One filepath per reviewer spreadsheet. Identical rows across spreadsheets are
# preserved as independent classifications.
CLASSIFIER_FILEPATHS = [
    Path(r"<YOUR PATH HERE>/S4.BlindAnomalyClassification_R1.xlsx"),
    Path(r"<YOUR PATH HERE>/S4.BlindAnomalyClassification_R2.xlsx"),
    Path(r"<YOUR PATH HERE>/S4.BlindAnomalyClassification_R3.xlsx"),
    Path(r"<YOUR PATH HERE>/S4.BlindAnomalyClassification_R4.xlsx")
]

# Optional RTK-guided Survey123 GeoPackage.
# The script will first try this configured Windows path, then fall back to a
# GeoPackage with the same filename in the current working directory, next to
# this script, or in /mnt/data when run in a notebook/sandbox.
INCLUDE_RTK_GUIDED_EVIDENCE = True
RTK_GPKG_PATH = Path(
    r"<YOUR PATH HERE>/S4.STADS_RTK_Anomalies_All_Merged_epsg4326.gpkg"
)
RTK_GPKG_FILENAME = "STADS_RTK_Anomalies_All_Merged_epsg4326.gpkg"
RTK_GPKG_LAYER = "Anomalies With Observations"
RTK_ANOMALIES_LAYER = "Anomalies"

# Optional RTK point metadata file. Use this if RTK sample_pai/a_or_b IDs do not
# exist in FILEPATH_FOR_LOCATIONS. The file should contain one row per RTK
# point, with columns for sample_pai, a_or_b, my_in_or_out, and stads_direction.
# It may be CSV, XLSX, or GPKG. Leave as None to fall back to pair-level RTK
# confirmation from compare_to only when no STADS lookup matches are found.
RTK_POINT_METADATA_PATH = None
RTK_POINT_METADATA_FILENAME = ""
RTK_POINT_METADATA_LAYER = None
RTK_POINT_METADATA_POINT_ID_COL = "sample_pai"
RTK_POINT_METADATA_PAIR_COL = "a_or_b"
RTK_POINT_METADATA_IN_OR_OUT_COL = "my_in_or_out"
RTK_POINT_METADATA_STADS_DIRECTION_COL = "stads_direction"

# If no RTK rows can be matched to a point-level STADS lookup, create one
# pair-level RTK confirmation record per sample_pai using compare_to only.
# This supports anomaly-presence summaries but cannot support direction agreement
# because the in-anomaly A/B point and STADS direction are unknown.
ALLOW_RTK_PAIR_LEVEL_FALLBACK_WITHOUT_STADS_LOOKUP = True
RTK_PAIR_LEVEL_LABEL = "PAIR_LEVEL"

# Columns present in the joined RTK/anomaly GeoPackage. A non-missing value in
# any of these columns indicates that the Survey123 point intersected a STADS
# anomaly polygon in the joined layer. For RTK-guided analyses, STADS anomaly
# direction and anomaly severity are evaluated only from z-scores. The STADS
# z-score convention used here is inverted relative to the intuitive crop-response
# sign: negative z-scores indicate positive crop deviance and positive z-scores
# indicate negative crop deviance. Residuals are intentionally not used for RTK
# anomaly direction or severity analyses.
RTK_JOINED_ANOMALY_ID_COLS = ["OBJECTID", "point_id", "path"]
RTK_JOINED_DIRECTION_NUMERIC_COL_CANDIDATES = ["zscore"]
RTK_JOINED_DIRECTION_SIGN_RULES = {
    "zscore": "inverse",
    "stads_zscore": "inverse",
}
RTK_JOINED_AREA_COL = "area"
RTK_JOINED_PERIMETER_COL = "perimeter"

# RTK / Survey123 field mapping from the uploaded GeoPackage.
# Survey123 truncates sample_pair_id to sample_pai in this GeoPackage.
# This is the point-pair identifier, not an individual point identifier.
RTK_POINT_ID_COL = "sample_pai"
RTK_PAIR_COL = "a_or_b"
RTK_CROP_COL = "crop"
RTK_CROP_OTHER_COL = "crop_other"
RTK_COMPARE_TO_COL = "compare_to"
RTK_OBSERVATION_COL = "crop_obser"
RTK_ISSUE_COLS = ["bad_things", "bad_thin_1"]
RTK_DATE_COL_CANDIDATES = ["CreationDa", "untitled_q", "EditDate"]
RTK_METHOD_LABEL = "RTK guided"
BLIND_METHOD_LABEL = "Blind paired"
RTK_CLASSIFIER_ID = "PL_RTK_field"

# Main output folder. A subdirectory is created so this does not overwrite the
# existing plotting outputs.
OUTPUT_DIR = Path(
    r"<YOUR PATH HERE>/ReviewerAssessment_Analysis"
)

# Core column names from the existing workflow.
POINT_ID_COL_NAME = "my_pair_id"
PAIR_ID_COL_NAME = "my_pair_a_or_b"
IN_OR_OUT_COL_NAME = "my_in_or_out"
ACTUAL_DIRECTION_COL_NAME = "stads_direction"
DIRECTION_INDICATED_COL_NAME = "Direction Indicated"

# Reviewer template columns. Extra spreadsheet columns are retained only if they
# are listed in OPTIONAL_CONTEXT_COLS or detected as evaluation-method columns.
EXPECTED_CLASSIFIER_COLS = [
    "Crop",
    "PointObservations",
    "PairedPointObservations",
    "Anomaly Indicated",
    "Anomaly Classification Confidence",
    "Classification Rationale",
    "Direction Indicated",
    "Direction Confidence",
    "Direction Classification Rationale",
    "row_id",
    "my_pair_id",
    "my_pair_a_or_b",
    "fn",
    "farmnumber",
    "fieldname",
    "sent_lon",
    "sent_lat",
]

# Optional context variables that are useful for modeling if present in either
# the classifier spreadsheets or the location/evaluation CSV. Edit these names to
# match your actual anomaly metric columns if available.
OPTIONAL_CONTEXT_COLS = [
    "evaluation_method",
    "Evaluation Method",
    "scouting_method",
    "Scouting Method",
    "dataset",
    "Dataset",
    "source",
    "Source",
    "anomaly_area_m2",
    "anomaly_area_ha",
    "area_m2",
    "area_ha",
    "stads_magnitude",
    "anomaly_magnitude",
    "magnitude",
    "ndvi_difference",
    "ndvi_diff",
    "mean_deviation",
    "max_deviation",
    "persistence",
    "days_detected",
    "image_date",
    "imagery_date",
    "date",
    "crop_stage",
    # RTK / Survey123 context columns
    "globalid",
    "CreationDa",
    "Creator",
    "EditDate",
    "Editor",
    "sample_pai",
    "a_or_b",
    "crop",
    "crop_other",
    "compare_to",
    "bad_things",
    "bad_thin_1",
    "crop_obser",
    "untitled_q",
    "photo_index_wide_photo_1_path",
    "photo_index_wide_photo_2_path",
    "photo_index_wide_photo_3_path",
    "photo_index_wide_photo_4_path",
    "X",
    "Y",
    # Joined STADS anomaly attributes from STADS_RTK_Anomalies_All_Merged_epsg4326.gpkg
    "fid_2",
    "OBJECTID",
    "anomaly_group",
    "point_id",
    "ndvi",
    "solus_pc1",
    "elevation",
    "aspect",
    "slope",
    "tpi",
    "tri",
    "lat",
    "long",
    "Predicted_NDVI",
    "zscore",
    "abs_zscore",
    "stads_zscore",
    "abs_stads_zscore",
    "CropName",
    "area",
    "perimeter",
    "centroid",
    "imagerydate",
    "window_id",
    "cum_gdd",
    "cum_precip",
    "layer",
    "path",
    "anomaly_area_m2",
    "anomaly_perimeter_m",
    "rtk_joined_to_anomaly",
    "rtk_joined_anomaly_id",
    "rtk_joined_anomaly_point_id",
    "rtk_lookup_source",
    "rtk_stads_direction_source_col",
    "rtk_stads_direction_sign_rule",
]

EVALUATION_METHOD_COL_CANDIDATES = [
    "evaluation_method_analysis",
    "evaluation_method",
    "Evaluation Method",
    "scouting_method",
    "Scouting Method",
    "dataset",
    "Dataset",
    "source",
    "Source",
]

# Consensus rule parameters.
MIN_RATERS_FOR_CONSENSUS = 2
CONFIDENCE_WEIGHTED_SUPPORT_THRESHOLD = 0.50

# Category orders used for plots and ordinal scoring.
ANOMALY_INDICATED_ORDER = [
    "Confirmed",
    "Possibly Confirmed",
    "Not Confirmed",
    "Indeterminate",
    "Excluded",
]

CONFIDENCE_ORDER = [
    "High",
    "Medium",
    "Low",
    "Indeterminate",
]

DIRECTION_ORDER = [
    "Positive",
    "Negative",
    "Same",
    "Indeterminate",
    "Excluded",
]

# Label scoring. These are intentionally explicit so the analysis is auditable.
ANOMALY_SUPPORT_SCORE = {
    "Confirmed": 1.0,
    "Possibly Confirmed": 0.5,
    "Not Confirmed": 0.0,
}

CONFIDENCE_WEIGHT = {
    "High": 1.0,
    "Medium": 2.0 / 3.0,
    "Low": 1.0 / 3.0,
    "Indeterminate": np.nan,
}

# Matplotlib export settings.
SAVE_FIGURES = True
FIGURE_DPI = 600


# =============================================================================
# GENERAL HELPERS
# =============================================================================


def clean_key_col(series: pd.Series) -> pd.Series:
    """Clean key columns while preserving missing values as pd.NA.

    Survey123 / GeoPackage exports often store integer IDs as floating point
    values such as 9906.0, while reviewer spreadsheets may store the same ID as
    9906. This helper normalizes integer-like float strings to their integer
    representation before constructing composite keys.
    """

    def _clean_one(value):
        if pd.isna(value):
            return pd.NA
        text = str(value).strip()
        if text in {"", "nan", "None", "<NA>", "NaN"}:
            return pd.NA
        try:
            as_float = float(text)
            if as_float.is_integer() and text.replace(".", "", 1).replace("-", "", 1).isdigit():
                return str(int(as_float))
        except Exception:
            pass
        return text

    return series.map(_clean_one)


def clean_text_col(series: pd.Series) -> pd.Series:
    """Clean free-text categorical columns."""
    return (
        series.astype(str)
        .str.strip()
        .replace({"nan": pd.NA, "None": pd.NA, "<NA>": pd.NA, "": pd.NA})
    )


def first_non_missing(values: Iterable) -> object:
    """Return the first non-missing value from an iterable, otherwise pd.NA."""
    for value in values:
        if pd.notna(value):
            return value
    return pd.NA


def mode_or_tie(values: Iterable) -> object:
    """Return the unique mode or 'Tie' if multiple values share the top count."""
    cleaned = pd.Series(list(values)).dropna()
    if cleaned.empty:
        return pd.NA
    counts = cleaned.value_counts()
    top = counts[counts == counts.max()].index.tolist()
    return top[0] if len(top) == 1 else "Tie"


def threshold_count(n_valid: int, threshold_fraction: float) -> int:
    """Number of ratings required to meet a consensus fraction."""
    if n_valid <= 0:
        return 0
    return int(math.ceil(threshold_fraction * n_valid))


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    """Return proportion, lower CI, upper CI using Wilson score interval."""
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half_width = z * math.sqrt((p * (1 - p) / n) + (z**2 / (4 * n**2))) / denom
    return (p, max(0.0, center - half_width), min(1.0, center + half_width))


def exact_binomial_p_value_greater_equal(k: int, n: int, p0: float = 0.5) -> float:
    """One-sided exact binomial p-value P[X >= k] under p0."""
    if n <= 0 or k < 0:
        return np.nan
    return float(
        sum(
            math.comb(n, x) * (p0**x) * ((1 - p0) ** (n - x))
            for x in range(k, n + 1)
        )
    )


def summarize_binary_rate(
    data: pd.DataFrame,
    success_col: str,
    group_cols: Optional[list[str]] = None,
    label: str = "rate",
) -> pd.DataFrame:
    """Summarize binary success rates with Wilson confidence intervals."""
    group_cols = group_cols or []
    rows = []

    if group_cols:
        grouped = data.groupby(group_cols, dropna=False)
    else:
        grouped = [((), data)]

    for group_key, group in grouped:
        success = group[success_col].dropna().astype(bool)
        n = int(len(success))
        k = int(success.sum())
        prop, lower, upper = wilson_ci(k, n)

        row = {
            "n": n,
            "successes": k,
            label: prop,
            f"{label}_percent": prop * 100 if pd.notna(prop) else np.nan,
            "ci_lower": lower,
            "ci_upper": upper,
            "ci_lower_percent": lower * 100 if pd.notna(lower) else np.nan,
            "ci_upper_percent": upper * 100 if pd.notna(upper) else np.nan,
        }

        if group_cols:
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            for col, value in zip(group_cols, group_key):
                row[col] = value

        rows.append(row)

    cols = group_cols + [
        "n",
        "successes",
        label,
        f"{label}_percent",
        "ci_lower",
        "ci_upper",
        "ci_lower_percent",
        "ci_upper_percent",
    ]
    return pd.DataFrame(rows)[cols]


# =============================================================================
# NORMALIZATION HELPERS
# =============================================================================


def normalize_anomaly_value(value) -> object:
    """Normalize anomaly-confirmation labels to a small set of categories."""
    if pd.isna(value):
        return pd.NA

    raw = str(value).strip()
    value_cf = raw.casefold()

    if value_cf in {"confirmed", "confirm", "yes", "y", "true"}:
        return "Confirmed"

    if value_cf in {
        "possibly confirmed",
        "possible confirmed",
        "possible",
        "possibly",
        "maybe",
        "partially confirmed",
        "partial",
    }:
        return "Possibly Confirmed"

    if value_cf in {
        "not confirmed",
        "not-confirmed",
        "no",
        "n",
        "false",
        "not",
        "unconfirmed",
    }:
        return "Not Confirmed"

    if value_cf in {"indeterminate", "unknown", "unclear", "na", "n/a", "nan", ""}:
        return "Indeterminate"

    if value_cf in {"excluded", "exclude", "bad data", "invalid"}:
        return "Excluded"

    return raw


def normalize_confidence_value(value) -> object:
    """Normalize confidence labels."""
    if pd.isna(value):
        return pd.NA

    raw = str(value).strip()
    value_cf = raw.casefold()

    if value_cf in {"high", "h"}:
        return "High"
    if value_cf in {"medium", "med", "m"}:
        return "Medium"
    if value_cf in {"low", "l"}:
        return "Low"
    if value_cf in {"indeterminate", "unknown", "unclear", "na", "n/a", "nan", ""}:
        return "Indeterminate"

    return raw


def normalize_direction_value(value) -> object:
    """Normalize direction labels to Positive / Negative / Same / Indeterminate."""
    if pd.isna(value):
        return pd.NA

    raw = str(value).strip()
    value_cf = raw.casefold()

    if value_cf in {"positive", "pos", "+", "higher", "better", "above", "high"}:
        return "Positive"

    if value_cf in {"negative", "neg", "-", "lower", "worse", "below", "low"}:
        return "Negative"

    if value_cf in {"same", "neutral", "no difference", "none", "similar"}:
        return "Same"

    if value_cf in {"indeterminate", "unknown", "unclear", "na", "n/a", "nan", ""}:
        return "Indeterminate"

    if value_cf in {"excluded", "exclude", "bad data", "invalid"}:
        return "Excluded"

    return raw


# =============================================================================
# DATA LOADING
# =============================================================================


def build_location_lookup() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load location/evaluation CSV and return full data plus in-anomaly lookup."""
    df_locations = pd.read_csv(FILEPATH_FOR_LOCATIONS)

    required_location_cols = [
        POINT_ID_COL_NAME,
        PAIR_ID_COL_NAME,
        IN_OR_OUT_COL_NAME,
        ACTUAL_DIRECTION_COL_NAME,
    ]
    missing_location_cols = [
        col for col in required_location_cols if col not in df_locations.columns
    ]
    if missing_location_cols:
        raise ValueError(
            f"The location file is missing required columns: {missing_location_cols}"
        )

    df_locations = df_locations.copy()
    for col in [POINT_ID_COL_NAME, PAIR_ID_COL_NAME, IN_OR_OUT_COL_NAME, ACTUAL_DIRECTION_COL_NAME]:
        df_locations[col] = clean_key_col(df_locations[col])

    df_locations["_composite_point_key"] = (
        df_locations[POINT_ID_COL_NAME].astype(str)
        + "||"
        + df_locations[PAIR_ID_COL_NAME].astype(str)
    )

    # Check for conflicting in/out labels.
    location_key_status_counts = (
        df_locations.dropna(subset=["_composite_point_key", IN_OR_OUT_COL_NAME])
        .groupby("_composite_point_key")[IN_OR_OUT_COL_NAME]
        .nunique()
    )
    conflicting_location_keys = location_key_status_counts[
        location_key_status_counts > 1
    ].index.tolist()
    if conflicting_location_keys:
        print(
            f"WARNING: {len(conflicting_location_keys)} composite keys have conflicting "
            f"{IN_OR_OUT_COL_NAME} values in the location file."
        )

    # Check for conflicting direction labels.
    location_key_direction_counts = (
        df_locations.dropna(subset=["_composite_point_key", ACTUAL_DIRECTION_COL_NAME])
        .groupby("_composite_point_key")[ACTUAL_DIRECTION_COL_NAME]
        .nunique()
    )
    conflicting_direction_keys = location_key_direction_counts[
        location_key_direction_counts > 1
    ].index.tolist()
    if conflicting_direction_keys:
        print(
            f"WARNING: {len(conflicting_direction_keys)} composite keys have conflicting "
            f"{ACTUAL_DIRECTION_COL_NAME} values in the location file."
        )

    context_cols = [
        col for col in OPTIONAL_CONTEXT_COLS if col in df_locations.columns
    ]
    lookup_cols = ["_composite_point_key", ACTUAL_DIRECTION_COL_NAME] + context_cols
    lookup_cols = list(dict.fromkeys(lookup_cols))

    in_anomaly_lookup = (
        df_locations[
            df_locations[IN_OR_OUT_COL_NAME].str.casefold() == "in"
        ][lookup_cols]
        .dropna(subset=["_composite_point_key"])
        .drop_duplicates(subset=["_composite_point_key"])
    )

    print(f"Number of unique in-anomaly point keys: {len(in_anomaly_lookup)}")
    return df_locations, in_anomaly_lookup


def load_reviewer_spreadsheets(in_anomaly_lookup: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load reviewer spreadsheets, filter to in-anomaly rows, and combine."""
    dfs = []
    filtered_out_dfs = []

    for path in CLASSIFIER_FILEPATHS:
        temp_df = pd.read_excel(path)

        missing_cols = [col for col in EXPECTED_CLASSIFIER_COLS if col not in temp_df.columns]
        extra_cols = [col for col in temp_df.columns if col not in EXPECTED_CLASSIFIER_COLS]

        if missing_cols:
            print(f"WARNING: {path.name} is missing columns: {missing_cols}")
        if extra_cols:
            print(f"WARNING: {path.name} has extra columns: {extra_cols}")

        required_classifier_key_cols = [POINT_ID_COL_NAME, PAIR_ID_COL_NAME]
        missing_classifier_key_cols = [
            col for col in required_classifier_key_cols if col not in temp_df.columns
        ]
        if missing_classifier_key_cols:
            raise ValueError(
                f"{path.name} is missing required key columns needed for in-anomaly filtering: "
                f"{missing_classifier_key_cols}"
            )

        temp_df = temp_df.copy()
        temp_df[POINT_ID_COL_NAME] = clean_key_col(temp_df[POINT_ID_COL_NAME])
        temp_df[PAIR_ID_COL_NAME] = clean_key_col(temp_df[PAIR_ID_COL_NAME])
        temp_df["_composite_point_key"] = (
            temp_df[POINT_ID_COL_NAME].astype(str)
            + "||"
            + temp_df[PAIR_ID_COL_NAME].astype(str)
        )

        rows_before_filter = len(temp_df)

        temp_df = temp_df.merge(
            in_anomaly_lookup.assign(_is_in_anomaly=True),
            on="_composite_point_key",
            how="left",
        )

        filtered_out_df = temp_df[temp_df["_is_in_anomaly"].isna()].copy()
        filtered_out_df["classifier_file"] = path.name
        filtered_out_df["classifier_file_path"] = str(path)

        temp_df = temp_df[temp_df["_is_in_anomaly"] == True].copy()

        if DIRECTION_INDICATED_COL_NAME not in temp_df.columns:
            raise ValueError(
                f"{path.name} is missing required direction classification column: "
                f"{DIRECTION_INDICATED_COL_NAME}"
            )

        rows_after_filter = len(temp_df)
        rows_removed = rows_before_filter - rows_after_filter
        print(
            f"{path.name}: kept {rows_after_filter} in-anomaly rows; "
            f"removed {rows_removed} out-of-anomaly or unmatched rows."
        )

        # Add classifier/source metadata so identical classifications remain traceable.
        temp_df["classifier_file"] = path.name
        temp_df["classifier_file_path"] = str(path)
        temp_df["classifier_id"] = (
            path.stem.replace("BlindAnomalyClassification_", "")
            .replace("AnomalyClassification_", "")
        )
        temp_df["evaluation_method"] = BLIND_METHOD_LABEL
        temp_df["evaluation_method_analysis"] = BLIND_METHOD_LABEL

        keep_cols = (
            EXPECTED_CLASSIFIER_COLS
            + OPTIONAL_CONTEXT_COLS
            + [
                ACTUAL_DIRECTION_COL_NAME,
                "_composite_point_key",
                "_is_in_anomaly",
                "classifier_file",
                "classifier_file_path",
                "classifier_id",
                "evaluation_method",
                "evaluation_method_analysis",
            ]
        )
        keep_cols = [col for col in dict.fromkeys(keep_cols) if col in temp_df.columns]
        dfs.append(temp_df[keep_cols].copy())

        if len(filtered_out_df) > 0:
            filtered_cols = [
                col
                for col in EXPECTED_CLASSIFIER_COLS
                + OPTIONAL_CONTEXT_COLS
                + ["classifier_file", "classifier_file_path", "_composite_point_key"]
                if col in filtered_out_df.columns
            ]
            filtered_out_dfs.append(filtered_out_df[filtered_cols].copy())

    if not dfs:
        raise ValueError("No classifier spreadsheets were loaded.")

    df = pd.concat(dfs, ignore_index=True)
    df_filtered_out = (
        pd.concat(filtered_out_dfs, ignore_index=True)
        if filtered_out_dfs
        else pd.DataFrame()
    )

    print(f"Combined in-anomaly classification rows: {len(df)}")
    print(f"Filtered-out rows: {len(df_filtered_out)}")
    print(f"Number of classifier files: {len(CLASSIFIER_FILEPATHS)}")
    print(df["classifier_file"].value_counts())

    return df, df_filtered_out



def resolve_rtk_gpkg_path() -> Optional[Path]:
    """Resolve the RTK-guided GeoPackage path from common locations."""
    candidates = [
        RTK_GPKG_PATH,
        Path.cwd() / RTK_GPKG_FILENAME,
    ]

    try:
        candidates.append(Path(__file__).resolve().parent / RTK_GPKG_FILENAME)
    except Exception:
        pass

    candidates.append(Path("/mnt/data") / RTK_GPKG_FILENAME)

    seen = set()
    for candidate in candidates:
        candidate = Path(candidate)
        if str(candidate) in seen:
            continue
        seen.add(str(candidate))
        if candidate.exists():
            return candidate

    return None


def derive_rtk_confirmation_from_inferred_direction(value) -> object:
    """Map an inferred RTK anomaly-point direction to an anomaly-presence label.

    Positive/Negative indicate a visible field difference at the in-anomaly point
    relative to its paired point and are treated as confirmed. Same is treated as
    not confirmed. Missing, tied, or ambiguous evidence is indeterminate.
    """
    direction = normalize_direction_value(value)
    if pd.isna(direction):
        return "Indeterminate"
    if direction in {"Positive", "Negative"}:
        return "Confirmed"
    if direction == "Same":
        return "Not Confirmed"
    return "Indeterminate"


def invert_pairwise_direction(direction) -> object:
    """Invert a current-point-vs-paired-point direction for the paired point.

    If the current point is Better/Positive than the paired point, then the
    paired point is Negative relative to the current point. Same remains Same.
    """
    direction = normalize_direction_value(direction)
    if pd.isna(direction):
        return pd.NA
    if direction == "Positive":
        return "Negative"
    if direction == "Negative":
        return "Positive"
    if direction == "Same":
        return "Same"
    return pd.NA


def opposite_pair_label(value) -> object:
    """Return the paired A/B label when available."""
    if pd.isna(value):
        return pd.NA
    text = str(value).strip().upper()
    if text == "A":
        return "B"
    if text == "B":
        return "A"
    return pd.NA


def is_in_anomaly_label(value) -> bool:
    """True when a location/evaluation label means the point is inside an anomaly."""
    if pd.isna(value):
        return False
    return str(value).strip().casefold() == "in"


def build_full_point_lookup(df_locations: pd.DataFrame) -> pd.DataFrame:
    """Build a lookup for all points, not only in-anomaly points.

    This is needed for RTK-guided scouting because compare_to may be entered on
    the out-of-anomaly point. In that case, the direction for the in-anomaly
    point is inferred by inversion.
    """
    context_cols = [col for col in OPTIONAL_CONTEXT_COLS if col in df_locations.columns]
    lookup_cols = [
        "_composite_point_key",
        POINT_ID_COL_NAME,
        PAIR_ID_COL_NAME,
        IN_OR_OUT_COL_NAME,
        ACTUAL_DIRECTION_COL_NAME,
    ] + context_cols
    lookup_cols = [col for col in dict.fromkeys(lookup_cols) if col in df_locations.columns]
    return (
        df_locations[lookup_cols]
        .dropna(subset=["_composite_point_key"])
        .drop_duplicates(subset=["_composite_point_key"])
        .copy()
    )


def get_numeric_direction_sign_rule(source_col: Optional[str]) -> str:
    """Return the STADS direction sign convention for a numeric anomaly metric.

    Current STADS RTK joined-anomaly outputs use an inverted z-score convention:
    negative z-scores indicate positive crop deviance, while positive z-scores
    indicate negative crop deviance. Residuals are intentionally not used for
    RTK anomaly direction or severity analyses.
    """
    if source_col is None or pd.isna(source_col):
        return "standard"
    source_col_cf = str(source_col).strip().casefold()
    return RTK_JOINED_DIRECTION_SIGN_RULES.get(source_col_cf, "standard")


def direction_from_numeric_anomaly_value(value, source_col: Optional[str] = None) -> object:
    """Map a signed STADS numeric anomaly value to Positive/Negative.

    For RTK-guided analyses, only zscore/stads_zscore should be passed here.
    The z-score sign convention is inverted: negative values are positive crop
    deviance and positive values are negative crop deviance.
    """
    if pd.isna(value):
        return pd.NA
    try:
        numeric_value = float(value)
    except Exception:
        return pd.NA

    sign_rule = get_numeric_direction_sign_rule(source_col)

    if numeric_value == 0:
        return pd.NA

    if sign_rule == "inverse":
        return "Positive" if numeric_value < 0 else "Negative"

    return "Positive" if numeric_value > 0 else "Negative"


def derive_joined_rtk_anomaly_fields(rtk: pd.DataFrame) -> pd.DataFrame:
    """Derive in/out and STADS-direction fields from the joined RTK/anomaly layer.

    The joined GeoPackage supplied for the RTK-guided analysis already contains
    Survey123 observation points spatially joined to STADS anomaly polygons. This
    helper converts those joined anomaly attributes into the same fields used by
    the rest of the analysis: my_in_or_out and stads_direction.
    """
    rtk = rtk.copy()
    id_cols = [col for col in RTK_JOINED_ANOMALY_ID_COLS if col in rtk.columns]
    if id_cols:
        joined_mask = rtk[id_cols].notna().any(axis=1)
    else:
        joined_mask = pd.Series(False, index=rtk.index)

    rtk["rtk_joined_to_anomaly"] = joined_mask
    rtk[IN_OR_OUT_COL_NAME] = np.where(joined_mask, "in", "out")

    # Use only zscore for RTK STADS direction. Do not fall back to residual,
    # because residuals are not used for RTK direction or severity evaluation.
    direction_source_col = None
    for col in RTK_JOINED_DIRECTION_NUMERIC_COL_CANDIDATES:
        if col in rtk.columns and pd.to_numeric(rtk[col], errors="coerce").notna().any():
            direction_source_col = col
            break

    if direction_source_col is not None:
        direction_source = pd.to_numeric(rtk[direction_source_col], errors="coerce")
        rtk[ACTUAL_DIRECTION_COL_NAME] = direction_source.apply(
            lambda value: direction_from_numeric_anomaly_value(
                value,
                source_col=direction_source_col,
            )
        )
        rtk["rtk_stads_direction_source_col"] = direction_source_col
        rtk["rtk_stads_direction_sign_rule"] = get_numeric_direction_sign_rule(direction_source_col)
    else:
        rtk[ACTUAL_DIRECTION_COL_NAME] = pd.NA
        rtk["rtk_stads_direction_source_col"] = pd.NA
        rtk["rtk_stads_direction_sign_rule"] = pd.NA

    if "zscore" in rtk.columns:
        rtk["stads_zscore"] = pd.to_numeric(rtk["zscore"], errors="coerce")
        rtk["abs_stads_zscore"] = rtk["stads_zscore"].abs()
        rtk["abs_zscore"] = rtk["stads_zscore"].abs()
    # Residuals are intentionally not converted into STADS severity fields here;
    # RTK anomaly direction and severity are evaluated only using z-scores.
    if RTK_JOINED_AREA_COL in rtk.columns:
        rtk["anomaly_area_m2"] = pd.to_numeric(rtk[RTK_JOINED_AREA_COL], errors="coerce")
    if RTK_JOINED_PERIMETER_COL in rtk.columns:
        rtk["anomaly_perimeter_m"] = pd.to_numeric(rtk[RTK_JOINED_PERIMETER_COL], errors="coerce")
    if "OBJECTID" in rtk.columns:
        rtk["rtk_joined_anomaly_id"] = rtk["OBJECTID"]
    if "point_id" in rtk.columns:
        rtk["rtk_joined_anomaly_point_id"] = rtk["point_id"]

    return rtk


def build_rtk_joined_point_lookup(rtk: pd.DataFrame) -> pd.DataFrame:
    """Build RTK point lookup directly from the joined RTK/anomaly GeoPackage.

    This is preferred over the blind location CSV when the RTK file already
    contains joined anomaly attributes. Duplicate point keys are retained when a
    point intersects multiple anomaly polygons, so ambiguous cases remain visible
    in the pairwise audit and collapse to Indeterminate when directions conflict.
    """
    required = ["_composite_point_key", POINT_ID_COL_NAME, PAIR_ID_COL_NAME, IN_OR_OUT_COL_NAME, ACTUAL_DIRECTION_COL_NAME]
    if any(col not in rtk.columns for col in required):
        return pd.DataFrame()

    if not rtk[IN_OR_OUT_COL_NAME].astype(str).str.casefold().eq("in").any():
        return pd.DataFrame()

    context_cols = [col for col in OPTIONAL_CONTEXT_COLS if col in rtk.columns]
    lookup_cols = required + context_cols
    lookup_cols = [col for col in dict.fromkeys(lookup_cols) if col in rtk.columns]
    out = rtk[lookup_cols].dropna(subset=["_composite_point_key"]).copy()
    out["rtk_lookup_source"] = "rtk_joined_anomaly_gpkg"
    return out


def make_rtk_base_dataframe(gdf: pd.DataFrame) -> pd.DataFrame:
    """Normalize the raw Survey123 GeoPackage rows before pairwise inference."""
    rtk = pd.DataFrame(gdf.drop(columns="geometry", errors="ignore")).copy()
    rtk["_rtk_raw_row_number"] = np.arange(len(rtk)) + 1

    required_cols = [RTK_POINT_ID_COL, RTK_PAIR_COL, RTK_COMPARE_TO_COL]
    missing_cols = [col for col in required_cols if col not in rtk.columns]
    if missing_cols:
        raise ValueError(
            f"The RTK GeoPackage is missing required columns: {missing_cols}. "
            "Update RTK_* column names in the CONFIG section."
        )

    rtk[POINT_ID_COL_NAME] = clean_key_col(rtk[RTK_POINT_ID_COL])
    rtk[PAIR_ID_COL_NAME] = clean_key_col(rtk[RTK_PAIR_COL])
    rtk["_composite_point_key"] = (
        rtk[POINT_ID_COL_NAME].astype(str) + "||" + rtk[PAIR_ID_COL_NAME].astype(str)
    )
    rtk["_paired_pair_label"] = rtk[PAIR_ID_COL_NAME].apply(opposite_pair_label)
    rtk["_paired_composite_point_key"] = (
        rtk[POINT_ID_COL_NAME].astype(str) + "||" + rtk["_paired_pair_label"].astype(str)
    )

    rtk["Crop"] = clean_text_col(rtk[RTK_CROP_COL]) if RTK_CROP_COL in rtk.columns else pd.NA
    if RTK_CROP_OTHER_COL in rtk.columns:
        crop_other = clean_text_col(rtk[RTK_CROP_OTHER_COL])
        other_mask = rtk["Crop"].astype(str).str.casefold().eq("other") & crop_other.notna()
        rtk.loc[other_mask, "Crop"] = crop_other[other_mask]

    rtk["PointObservations"] = (
        rtk[RTK_OBSERVATION_COL] if RTK_OBSERVATION_COL in rtk.columns else pd.NA
    )
    rtk["PairedPointObservations"] = pd.NA

    if "geometry" in gdf:
        try:
            rtk["sent_lon"] = gdf.geometry.x
            rtk["sent_lat"] = gdf.geometry.y
        except Exception:
            pass
    if "sent_lon" not in rtk.columns and "X" in rtk.columns:
        rtk["sent_lon"] = rtk["X"]
    if "sent_lat" not in rtk.columns and "Y" in rtk.columns:
        rtk["sent_lat"] = rtk["Y"]

    for date_col in RTK_DATE_COL_CANDIDATES:
        if date_col in rtk.columns:
            rtk["image_date"] = rtk[date_col]
            break

    # If this is the joined RTK/anomaly layer, use its anomaly attributes to
    # derive in/out status, STADS direction, and anomaly metric covariates.
    rtk = derive_joined_rtk_anomaly_fields(rtk)

    return rtk


def build_rtk_pairwise_evidence(
    rtk: pd.DataFrame,
    full_point_lookup: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Infer evidence for in-anomaly RTK points from pair-relative compare_to.

    Survey123 sample_pai is the point-pair identifier and a_or_b identifies the
    A or B point within that same pair. compare_to is entered on a visited/current
    row and describes that current A/B point relative to the other A/B point with
    the same sample_pai value. Therefore:

    - If compare_to is entered on the in-anomaly point, Better/Worse/Same is used
      directly for that in-anomaly point.
    - If compare_to is entered on the paired out-of-anomaly point, the direction
      is inverted before assigning evidence to the in-anomaly point.
    - If both A and B have compare_to values, both become pairwise evidence rows
      and are later collapsed to one RTK evidence row per in-anomaly point.
    """
    point_lookup = full_point_lookup.copy()

    current_lookup = point_lookup.add_prefix("current_")
    paired_lookup = point_lookup.add_prefix("paired_")

    enriched = rtk.merge(
        current_lookup,
        left_on="_composite_point_key",
        right_on="current__composite_point_key",
        how="left",
    ).merge(
        paired_lookup,
        left_on="_paired_composite_point_key",
        right_on="paired__composite_point_key",
        how="left",
    )

    evidence_rows = []
    no_target_rows = []

    for _, row in enriched.iterrows():
        raw_compare_to = row.get(RTK_COMPARE_TO_COL, pd.NA)
        current_vs_pair_direction = normalize_direction_value(raw_compare_to)
        current_key = row.get("_composite_point_key", pd.NA)
        paired_key = row.get("_paired_composite_point_key", pd.NA)
        current_in = is_in_anomaly_label(row.get(f"current_{IN_OR_OUT_COL_NAME}", pd.NA))
        paired_in = is_in_anomaly_label(row.get(f"paired_{IN_OR_OUT_COL_NAME}", pd.NA))

        candidate_targets = []
        if current_in:
            candidate_targets.append(
                {
                    "target_composite_point_key": current_key,
                    "target_role": "current_row_is_in_anomaly",
                    "inferred_anomaly_direction": current_vs_pair_direction,
                    "direction_was_inverted": False,
                    "target_stads_direction": row.get(f"current_{ACTUAL_DIRECTION_COL_NAME}", pd.NA),
                    "target_pair_label": row.get(PAIR_ID_COL_NAME, pd.NA),
                }
            )
        if paired_in:
            candidate_targets.append(
                {
                    "target_composite_point_key": paired_key,
                    "target_role": "paired_row_is_in_anomaly_direction_inverted",
                    "inferred_anomaly_direction": invert_pairwise_direction(current_vs_pair_direction),
                    "direction_was_inverted": True,
                    "target_stads_direction": row.get(f"paired_{ACTUAL_DIRECTION_COL_NAME}", pd.NA),
                    "target_pair_label": row.get("_paired_pair_label", pd.NA),
                }
            )

        if not candidate_targets:
            no_target = row.to_dict()
            no_target["rtk_pairwise_status"] = "no_in_anomaly_point_for_current_or_paired_key"
            no_target_rows.append(no_target)
            continue

        for candidate in candidate_targets:
            target_key = candidate["target_composite_point_key"]
            if pd.isna(target_key) or str(target_key).endswith("||<NA>") or str(target_key).endswith("||nan"):
                no_target = row.to_dict()
                no_target["rtk_pairwise_status"] = "invalid_target_key"
                no_target_rows.append(no_target)
                continue

            point_id, pair_label = (str(target_key).split("||", 1) + [pd.NA])[:2]
            inferred_direction = candidate["inferred_anomaly_direction"]
            if pd.isna(inferred_direction) or inferred_direction not in {"Positive", "Negative", "Same"}:
                pairwise_status = "no_evaluable_compare_to"
            else:
                pairwise_status = "evaluable_pairwise_compare_to"

            evidence = {
                "target_composite_point_key": target_key,
                POINT_ID_COL_NAME: point_id,
                PAIR_ID_COL_NAME: pair_label,
                "source_composite_point_key": current_key,
                "source_pair_label": row.get(PAIR_ID_COL_NAME, pd.NA),
                "source_is_target_point": current_key == target_key,
                "rtk_target_role": candidate["target_role"],
                "direction_was_inverted": candidate["direction_was_inverted"],
                "rtk_pairwise_status": pairwise_status,
                "rtk_raw_compare_to": raw_compare_to,
                "current_vs_pair_direction_normalized": current_vs_pair_direction,
                "inferred_anomaly_direction": inferred_direction,
                "inferred_anomaly_confirmation": derive_rtk_confirmation_from_inferred_direction(inferred_direction),
                ACTUAL_DIRECTION_COL_NAME: candidate["target_stads_direction"],
                "Crop": row.get("Crop", pd.NA),
                "PointObservations": row.get("PointObservations", pd.NA),
                "PairedPointObservations": row.get("PairedPointObservations", pd.NA),
                "sent_lon": row.get("sent_lon", pd.NA),
                "sent_lat": row.get("sent_lat", pd.NA),
                "image_date": row.get("image_date", pd.NA),
                "_rtk_raw_row_number": row.get("_rtk_raw_row_number", pd.NA),
            }

            for col in OPTIONAL_CONTEXT_COLS:
                if col in row.index and col not in evidence:
                    evidence[col] = row.get(col, pd.NA)

            evidence_rows.append(evidence)

    evidence_df = pd.DataFrame(evidence_rows)
    no_target_df = pd.DataFrame(no_target_rows)
    return evidence_df, no_target_df


def collapse_rtk_pairwise_evidence_to_points(pairwise: pd.DataFrame, gpkg_path: Path) -> pd.DataFrame:
    """Collapse possibly multiple RTK pairwise rows to one row per anomaly point."""
    if pairwise.empty:
        return pd.DataFrame()

    rows = []
    for target_key, group in pairwise.groupby("target_composite_point_key", dropna=False):
        valid_dirs = group[
            group["inferred_anomaly_direction"].isin(["Positive", "Negative", "Same"])
        ]["inferred_anomaly_direction"]
        direction_mode = mode_or_tie(valid_dirs)
        if pd.isna(direction_mode):
            direction_for_analysis = "Indeterminate"
        elif direction_mode == "Tie":
            direction_for_analysis = "Indeterminate"
        else:
            direction_for_analysis = direction_mode

        confirmation_labels = group[
            group["inferred_anomaly_confirmation"].isin(["Confirmed", "Not Confirmed"])
        ]["inferred_anomaly_confirmation"]
        confirmation_mode = mode_or_tie(confirmation_labels)
        if pd.isna(confirmation_mode):
            confirmation_for_analysis = "Indeterminate"
        elif confirmation_mode == "Tie":
            confirmation_for_analysis = "Indeterminate"
        else:
            confirmation_for_analysis = confirmation_mode

        n_valid = int(valid_dirs.notna().sum())
        n_total = int(len(group))
        n_distinct_valid_dirs = int(valid_dirs.nunique()) if n_valid else 0
        n_inverted = int(group["direction_was_inverted"].fillna(False).astype(bool).sum())

        if n_valid == 0:
            confidence = "Indeterminate"
        elif direction_for_analysis == "Indeterminate" or n_distinct_valid_dirs > 1:
            confidence = "Low"
        else:
            confidence = "High"

        point_id, pair_label = (str(target_key).split("||", 1) + [pd.NA])[:2]

        role_counts = group["rtk_target_role"].dropna().astype(str).value_counts().to_dict()
        compare_values = "; ".join(
            group["rtk_raw_compare_to"].dropna().astype(str).unique().tolist()
        )
        inferred_values = "; ".join(
            group["inferred_anomaly_direction"].dropna().astype(str).unique().tolist()
        )
        source_keys = "; ".join(
            group["source_composite_point_key"].dropna().astype(str).unique().tolist()
        )

        stads_direction_values = group[ACTUAL_DIRECTION_COL_NAME].dropna().astype(str) if ACTUAL_DIRECTION_COL_NAME in group.columns else pd.Series(dtype=object)
        stads_direction_mode = mode_or_tie(stads_direction_values)
        if pd.isna(stads_direction_mode) or stads_direction_mode == "Tie":
            stads_direction_for_analysis = pd.NA
        else:
            stads_direction_for_analysis = stads_direction_mode

        row = {
            "_composite_point_key": target_key,
            POINT_ID_COL_NAME: point_id,
            PAIR_ID_COL_NAME: pair_label,
            "Crop": first_non_missing(group["Crop"]),
            "PointObservations": "\n---\n".join(
                group["PointObservations"].dropna().astype(str).unique().tolist()
            ),
            "PairedPointObservations": pd.NA,
            "Anomaly Indicated": confirmation_for_analysis,
            "Anomaly Classification Confidence": confidence,
            "Classification Rationale": (
                "RTK-guided field assessment inferred from pair-relative Survey123 compare_to. "
                f"Evidence rows: {n_total}; evaluable compare_to rows: {n_valid}; "
                f"raw compare_to values: {compare_values or 'none'}; "
                f"inferred anomaly directions: {inferred_values or 'none'}. "
                "When compare_to was entered on the out-of-anomaly paired point, "
                "the direction was inverted before assigning evidence to the in-anomaly point."
            ),
            "Direction Indicated": direction_for_analysis,
            "Direction Confidence": confidence,
            "Direction Classification Rationale": (
                "Direction is the modal inferred direction for the in-anomaly point after "
                "accounting for whether compare_to was entered on the anomaly point or on "
                "the paired comparison point. Conflicting inferred directions are marked "
                "Indeterminate with Low confidence."
            ),
            ACTUAL_DIRECTION_COL_NAME: stads_direction_for_analysis,
            "rtk_distinct_stads_directions": int(stads_direction_values.nunique()) if len(stads_direction_values) else 0,
            "rtk_stads_direction_values": "; ".join(stads_direction_values.unique().tolist()) if len(stads_direction_values) else "",
            "rtk_joined_anomaly_ids": "; ".join(group["rtk_joined_anomaly_id"].dropna().astype(str).unique().tolist()) if "rtk_joined_anomaly_id" in group.columns else "",
            "sent_lon": first_non_missing(group["sent_lon"]) if "sent_lon" in group.columns else pd.NA,
            "sent_lat": first_non_missing(group["sent_lat"]) if "sent_lat" in group.columns else pd.NA,
            "classifier_file": gpkg_path.name,
            "classifier_file_path": str(gpkg_path),
            "classifier_id": RTK_CLASSIFIER_ID,
            "evaluation_method": RTK_METHOD_LABEL,
            "evaluation_method_analysis": RTK_METHOD_LABEL,
            "rtk_pairwise_evidence_rows": n_total,
            "rtk_pairwise_evaluable_rows": n_valid,
            "rtk_pairwise_distinct_valid_directions": n_distinct_valid_dirs,
            "rtk_pairwise_inverted_rows": n_inverted,
            "rtk_pairwise_role_counts": str(role_counts),
            "rtk_pairwise_source_keys": source_keys,
        }

        for col in OPTIONAL_CONTEXT_COLS:
            if col in group.columns and col not in row:
                row[col] = first_non_missing(group[col])

        rows.append(row)

    return pd.DataFrame(rows)




def resolve_optional_rtk_point_metadata_path() -> Optional[Path]:
    """Resolve optional RTK point metadata path from common locations."""
    if RTK_POINT_METADATA_PATH is None and not RTK_POINT_METADATA_FILENAME:
        return None

    candidates = []
    if RTK_POINT_METADATA_PATH is not None:
        candidates.append(Path(RTK_POINT_METADATA_PATH))
    if RTK_POINT_METADATA_FILENAME:
        candidates.append(Path.cwd() / RTK_POINT_METADATA_FILENAME)
        try:
            candidates.append(Path(__file__).resolve().parent / RTK_POINT_METADATA_FILENAME)
        except Exception:
            pass
        candidates.append(Path("/mnt/data") / RTK_POINT_METADATA_FILENAME)

    seen = set()
    for candidate in candidates:
        candidate = Path(candidate)
        if str(candidate) in seen:
            continue
        seen.add(str(candidate))
        if candidate.exists():
            return candidate
    return None


def load_optional_rtk_point_metadata() -> pd.DataFrame:
    """Load optional point-level RTK metadata with in/out and STADS direction.

    This is needed when RTK sample_pai/a_or_b values are not present in the
    blind matched-point location CSV. The metadata file should include the RTK
    point-pair identifier, the A/B label, an in/out label, and STADS direction.
    """
    metadata_path = resolve_optional_rtk_point_metadata_path()
    if metadata_path is None:
        return pd.DataFrame()

    suffix = metadata_path.suffix.casefold()
    if suffix == ".csv":
        meta = pd.read_csv(metadata_path)
    elif suffix in {".xlsx", ".xls"}:
        meta = pd.read_excel(metadata_path)
    elif suffix in {".gpkg", ".geojson", ".json", ".shp"}:
        try:
            import geopandas as gpd
        except Exception as exc:
            raise ImportError(
                "geopandas is required to read RTK point metadata spatial files."
            ) from exc
        if RTK_POINT_METADATA_LAYER:
            meta = gpd.read_file(metadata_path, layer=RTK_POINT_METADATA_LAYER)
        else:
            meta = gpd.read_file(metadata_path)
        meta = pd.DataFrame(meta.drop(columns="geometry", errors="ignore"))
    else:
        raise ValueError(
            f"Unsupported RTK point metadata file extension: {metadata_path.suffix}. "
            "Use CSV, XLSX, GPKG, GeoJSON, or SHP."
        )

    required = [
        RTK_POINT_METADATA_POINT_ID_COL,
        RTK_POINT_METADATA_PAIR_COL,
        RTK_POINT_METADATA_IN_OR_OUT_COL,
        RTK_POINT_METADATA_STADS_DIRECTION_COL,
    ]
    missing = [col for col in required if col not in meta.columns]
    if missing:
        raise ValueError(
            f"RTK point metadata file {metadata_path} is missing required columns: {missing}"
        )

    out = pd.DataFrame()
    out[POINT_ID_COL_NAME] = clean_key_col(meta[RTK_POINT_METADATA_POINT_ID_COL])
    out[PAIR_ID_COL_NAME] = clean_key_col(meta[RTK_POINT_METADATA_PAIR_COL])
    out[IN_OR_OUT_COL_NAME] = clean_key_col(meta[RTK_POINT_METADATA_IN_OR_OUT_COL])
    out[ACTUAL_DIRECTION_COL_NAME] = clean_key_col(meta[RTK_POINT_METADATA_STADS_DIRECTION_COL])
    out["_composite_point_key"] = (
        out[POINT_ID_COL_NAME].astype(str) + "||" + out[PAIR_ID_COL_NAME].astype(str)
    )

    for col in OPTIONAL_CONTEXT_COLS:
        if col in meta.columns and col not in out.columns:
            out[col] = meta[col]

    out = out.dropna(subset=["_composite_point_key"]).drop_duplicates(
        subset=["_composite_point_key"]
    )
    print(
        f"Loaded optional RTK point metadata from {metadata_path.name}: "
        f"{len(out)} unique point keys."
    )
    return out


def build_rtk_point_lookup(
    df_locations: pd.DataFrame,
    rtk_raw: Optional[pd.DataFrame] = None,
) -> tuple[pd.DataFrame, str]:
    """Return the best available point lookup for RTK in/out and STADS direction.

    Priority order:
    1. Joined RTK/anomaly GeoPackage attributes, when present.
    2. Optional external RTK metadata file.
    3. The blind matched-point location file as a final fallback.
    """
    if rtk_raw is not None:
        joined_lookup = build_rtk_joined_point_lookup(rtk_raw)
        if not joined_lookup.empty:
            return joined_lookup, "rtk_joined_anomaly_gpkg"

    rtk_metadata = load_optional_rtk_point_metadata()
    if not rtk_metadata.empty:
        return rtk_metadata, "optional_rtk_point_metadata"
    return build_full_point_lookup(df_locations), "main_location_file"


def summarize_rtk_key_overlap(rtk_raw: pd.DataFrame, point_lookup: pd.DataFrame) -> pd.DataFrame:
    """Summarize whether RTK current/paired keys occur in the STADS point lookup."""
    current_keys = set(rtk_raw["_composite_point_key"].dropna().astype(str))
    paired_keys = set(rtk_raw["_paired_composite_point_key"].dropna().astype(str))
    lookup_keys = set(point_lookup["_composite_point_key"].dropna().astype(str)) if not point_lookup.empty else set()
    return pd.DataFrame(
        [
            {
                "rtk_raw_rows": len(rtk_raw),
                "rtk_unique_current_keys": len(current_keys),
                "rtk_unique_paired_keys": len(paired_keys),
                "lookup_unique_keys": len(lookup_keys),
                "current_keys_in_lookup": len(current_keys & lookup_keys),
                "paired_keys_in_lookup": len(paired_keys & lookup_keys),
                "any_keys_in_lookup": len((current_keys | paired_keys) & lookup_keys),
            }
        ]
    )


def build_rtk_pair_level_fallback_evidence(
    rtk_raw: pd.DataFrame,
    gpkg_path: Path,
) -> pd.DataFrame:
    """Create pair-level RTK confirmation records when STADS point lookup is absent.

    This fallback uses compare_to only to answer whether the RTK-guided pair had
    an observed field difference. It does not assign an in-anomaly A/B point and
    does not evaluate STADS direction. Use point-level RTK metadata if you need
    RTK direction agreement with STADS.
    """
    if rtk_raw.empty:
        return pd.DataFrame()

    rows = []
    for pair_id, group in rtk_raw.groupby(POINT_ID_COL_NAME, dropna=False):
        valid_dirs = group[RTK_COMPARE_TO_COL].apply(normalize_direction_value)
        valid_dirs = valid_dirs[valid_dirs.isin(["Positive", "Negative", "Same"])]
        n_valid = int(valid_dirs.notna().sum())
        n_total = int(len(group))
        n_difference = int(valid_dirs.isin(["Positive", "Negative"]).sum())
        n_same = int((valid_dirs == "Same").sum())

        if n_valid == 0:
            anomaly_label = "Indeterminate"
            confidence = "Indeterminate"
        elif n_difference > 0:
            anomaly_label = "Confirmed"
            confidence = "High" if n_valid >= 1 else "Indeterminate"
        elif n_same == n_valid:
            anomaly_label = "Not Confirmed"
            confidence = "High" if n_valid >= 1 else "Indeterminate"
        else:
            anomaly_label = "Indeterminate"
            confidence = "Low"

        crop = first_non_missing(group["Crop"]) if "Crop" in group.columns else pd.NA
        obs_values = []
        if "PointObservations" in group.columns:
            obs_values = group["PointObservations"].dropna().astype(str).unique().tolist()

        raw_compare_values = "; ".join(
            group[RTK_COMPARE_TO_COL].dropna().astype(str).unique().tolist()
        ) if RTK_COMPARE_TO_COL in group.columns else ""
        normalized_compare_values = "; ".join(valid_dirs.dropna().astype(str).unique().tolist())

        row = {
            "_composite_point_key": f"{pair_id}||{RTK_PAIR_LEVEL_LABEL}",
            POINT_ID_COL_NAME: pair_id,
            PAIR_ID_COL_NAME: RTK_PAIR_LEVEL_LABEL,
            "Crop": crop,
            "PointObservations": "\n---\n".join(obs_values),
            "PairedPointObservations": pd.NA,
            "Anomaly Indicated": anomaly_label,
            "Anomaly Classification Confidence": confidence,
            "Classification Rationale": (
                "RTK pair-level fallback: no RTK sample_pai/a_or_b keys matched a "
                "point-level STADS lookup with in/out and direction labels. This row "
                "therefore summarizes whether the pair showed an observed difference "
                "from compare_to only. It supports anomaly-presence summaries, not "
                "STADS direction agreement."
            ),
            "Direction Indicated": "Indeterminate",
            "Direction Confidence": "Indeterminate",
            "Direction Classification Rationale": (
                "Direction agreement unavailable in RTK pair-level fallback because "
                "the in-anomaly A/B point and STADS direction were not available."
            ),
            ACTUAL_DIRECTION_COL_NAME: pd.NA,
            "sent_lon": pd.to_numeric(group["sent_lon"], errors="coerce").mean() if "sent_lon" in group.columns else pd.NA,
            "sent_lat": pd.to_numeric(group["sent_lat"], errors="coerce").mean() if "sent_lat" in group.columns else pd.NA,
            "classifier_file": gpkg_path.name,
            "classifier_file_path": str(gpkg_path),
            "classifier_id": RTK_CLASSIFIER_ID,
            "evaluation_method": RTK_METHOD_LABEL,
            "evaluation_method_analysis": RTK_METHOD_LABEL,
            "rtk_pairwise_evidence_rows": n_total,
            "rtk_pairwise_evaluable_rows": n_valid,
            "rtk_pairwise_distinct_valid_directions": int(valid_dirs.nunique()) if n_valid else 0,
            "rtk_pairwise_inverted_rows": 0,
            "rtk_pairwise_role_counts": "pair_level_fallback_without_stads_lookup",
            "rtk_pairwise_source_keys": "; ".join(group["_composite_point_key"].dropna().astype(str).unique().tolist()),
            "rtk_pair_level_fallback": True,
            "rtk_pair_level_valid_compare_to_rows": n_valid,
            "rtk_pair_level_difference_rows": n_difference,
            "rtk_pair_level_same_rows": n_same,
            "rtk_raw_compare_to_values": raw_compare_values,
            "rtk_normalized_compare_to_values": normalized_compare_values,
            "_is_in_anomaly": True,
        }

        for col in OPTIONAL_CONTEXT_COLS:
            if col in group.columns and col not in row:
                row[col] = first_non_missing(group[col])

        rows.append(row)

    return pd.DataFrame(rows)
def make_duplicate_rtk_analysis_keys_anomaly_specific(rtk_in: pd.DataFrame) -> pd.DataFrame:
    """Make RTK analysis keys anomaly-specific when one point intersects multiple anomalies.

    The Survey123 point key (sample_pai||a_or_b) can appear more than once in the
    joined RTK/anomaly layer if the same observation point intersects multiple
    STADS anomaly polygons. For RTK analysis, those rows should remain separate
    anomaly-specific evidence records so that direction agreement is evaluated
    against each joined anomaly's STADS direction instead of being collapsed into
    one ambiguous point key.
    """
    if rtk_in.empty or "_composite_point_key" not in rtk_in.columns:
        return rtk_in

    out = rtk_in.copy()
    dup_mask = out["_composite_point_key"].duplicated(keep=False)
    if not dup_mask.any():
        return out

    id_col = None
    for candidate in ["rtk_joined_anomaly_id", "OBJECTID", "point_id", "rtk_joined_anomaly_point_id"]:
        if candidate in out.columns and out.loc[dup_mask, candidate].notna().any():
            id_col = candidate
            break

    if id_col is None:
        out.loc[dup_mask, "_composite_point_key"] = (
            out.loc[dup_mask, "_composite_point_key"].astype(str)
            + "||rtk_duplicate_"
            + out.loc[dup_mask].groupby("_composite_point_key").cumcount().astype(str)
        )
    else:
        suffix = out.loc[dup_mask, id_col].astype(str).str.replace(r"\\.0$", "", regex=True)
        out.loc[dup_mask, "_composite_point_key"] = (
            out.loc[dup_mask, "_composite_point_key"].astype(str)
            + "||anom_"
            + suffix
        )

    # If the anomaly-id suffix still leaves duplicated keys, append a stable
    # within-key counter so downstream consensus does not collapse separate
    # joined-anomaly evidence rows.
    still_dup_mask = out["_composite_point_key"].duplicated(keep=False)
    if still_dup_mask.any():
        out.loc[still_dup_mask, "_composite_point_key"] = (
            out.loc[still_dup_mask, "_composite_point_key"].astype(str)
            + "||dup_"
            + out.loc[still_dup_mask].groupby("_composite_point_key").cumcount().astype(str)
        )

    out.loc[dup_mask, "rtk_analysis_key_was_made_anomaly_specific"] = True
    out.loc[~dup_mask, "rtk_analysis_key_was_made_anomaly_specific"] = False
    return out


def load_rtk_guided_gpkg(
    df_locations: pd.DataFrame,
    in_anomaly_lookup: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load RTK-guided Survey123 observations and infer anomaly evidence.

    The preferred RTK path is point-level inference using a lookup that contains
    RTK sample_pai/a_or_b, in/out status, and STADS direction. If the RTK IDs do
    not match that lookup, the function can fall back to pair-level evidence from
    compare_to only. Pair-level fallback supports confirmation-rate summaries but
    cannot support STADS direction agreement.
    """
    if not INCLUDE_RTK_GUIDED_EVIDENCE:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    gpkg_path = resolve_rtk_gpkg_path()
    if gpkg_path is None:
        print(
            "WARNING: RTK-guided GeoPackage was not found. "
            f"Edit RTK_GPKG_PATH or place {RTK_GPKG_FILENAME} next to this script."
        )
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    try:
        import geopandas as gpd
    except Exception as exc:
        print(
            "WARNING: geopandas is required to read the RTK-guided GeoPackage. "
            f"RTK evidence skipped. Import error: {repr(exc)}"
        )
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    try:
        gdf = gpd.read_file(gpkg_path, layer=RTK_GPKG_LAYER)
    except Exception:
        # Fall back to default/first layer if the layer name changes.
        gdf = gpd.read_file(gpkg_path)

    rtk_raw = make_rtk_base_dataframe(gdf)
    rtk_point_lookup, lookup_source = build_rtk_point_lookup(df_locations, rtk_raw=rtk_raw)
    key_overlap_summary = summarize_rtk_key_overlap(rtk_raw, rtk_point_lookup)
    print(
        "RTK key-overlap diagnostic: "
        + key_overlap_summary.to_dict(orient="records")[0].__repr__()
        + f"; lookup_source={lookup_source}"
    )

    pairwise_evidence, no_target_evidence = build_rtk_pairwise_evidence(
        rtk_raw,
        rtk_point_lookup,
    )
    rtk_point_evidence = collapse_rtk_pairwise_evidence_to_points(
        pairwise_evidence,
        gpkg_path,
    )

    # Keep only target points that are inside anomalies according to the lookup.
    rows_before_filter = len(rtk_point_evidence)
    if not rtk_point_evidence.empty:
        if lookup_source in {"optional_rtk_point_metadata", "rtk_joined_anomaly_gpkg"}:
            rtk_in_anomaly_lookup = rtk_point_lookup[
                rtk_point_lookup[IN_OR_OUT_COL_NAME].str.casefold() == "in"
            ].copy()
        else:
            rtk_in_anomaly_lookup = in_anomaly_lookup.copy()

        rtk_point_evidence = rtk_point_evidence.merge(
            rtk_in_anomaly_lookup.assign(_is_in_anomaly=True),
            on="_composite_point_key",
            how="left",
            suffixes=("", "_lookup"),
        )
        # If the merge introduced lookup copies of context columns, use them only
        # to fill missing values and then remove the duplicated columns.
        for col in list(rtk_point_evidence.columns):
            if col.endswith("_lookup"):
                base_col = col[:-7]
                if base_col in rtk_point_evidence.columns:
                    # For joined-anomaly attributes, prefer the lookup value so
                    # duplicated point keys created by multiple intersecting
                    # anomaly polygons keep the correct per-anomaly ID, sign,
                    # and metrics after the merge. For general context, keep
                    # existing values and use lookup only to fill gaps.
                    prefer_lookup_cols = {
                        ACTUAL_DIRECTION_COL_NAME,
                        "stads_direction_normalized",
                        "rtk_joined_anomaly_id",
                        "rtk_joined_anomaly_point_id",
                        "OBJECTID",
                        "point_id",
                        "anomaly_group",
                                            "zscore",
                        "stads_zscore",
                        "abs_stads_zscore",
                                                                "anomaly_area_m2",
                        "anomaly_perimeter_m",
                        "area",
                        "perimeter",
                        "CropName",
                        "imagerydate",
                        "layer",
                        "path",
                    }
                    if base_col in prefer_lookup_cols:
                        rtk_point_evidence[base_col] = rtk_point_evidence[col].combine_first(
                            rtk_point_evidence[base_col]
                        )
                    else:
                        rtk_point_evidence[base_col] = rtk_point_evidence[base_col].combine_first(
                            rtk_point_evidence[col]
                        )
                else:
                    rtk_point_evidence[base_col] = rtk_point_evidence[col]
                rtk_point_evidence = rtk_point_evidence.drop(columns=[col])

        rtk_filtered_out = rtk_point_evidence[rtk_point_evidence["_is_in_anomaly"].isna()].copy()
        rtk_in = rtk_point_evidence[rtk_point_evidence["_is_in_anomaly"] == True].copy()
    else:
        rtk_filtered_out = pd.DataFrame()
        rtk_in = pd.DataFrame()

    if not no_target_evidence.empty:
        no_target_evidence["_is_in_anomaly"] = pd.NA
        no_target_evidence["rtk_lookup_source"] = lookup_source
        rtk_filtered_out = pd.concat([rtk_filtered_out, no_target_evidence], ignore_index=True, sort=False)

    if rtk_in.empty and ALLOW_RTK_PAIR_LEVEL_FALLBACK_WITHOUT_STADS_LOOKUP:
        print(
            "WARNING: No RTK sample_pai/a_or_b keys matched a point-level STADS "
            "lookup with in/out labels. Creating RTK pair-level fallback evidence "
            "from compare_to only. This will support RTK confirmation-rate summaries "
            "but RTK direction agreement with STADS will remain unavailable."
        )
        rtk_in = build_rtk_pair_level_fallback_evidence(rtk_raw, gpkg_path)
        if not rtk_in.empty:
            rtk_in["rtk_lookup_source"] = "pair_level_fallback_without_stads_lookup"
        # Treat raw unmatched rows as an audit trail, but do not count them as
        # filtered-out analysis rows after fallback.
        if rtk_filtered_out.empty:
            rtk_filtered_out = rtk_raw.copy()
            rtk_filtered_out["rtk_lookup_source"] = lookup_source
            rtk_filtered_out["rtk_pairwise_status"] = "unmatched_before_pair_level_fallback"

    if not rtk_in.empty:
        rtk_in = make_duplicate_rtk_analysis_keys_anomaly_specific(rtk_in)

    print(
        f"RTK-guided GeoPackage {gpkg_path.name}: inferred {len(rtk_in)} analysis "
        f"evidence rows from {len(rtk_raw)} raw Survey123 rows; "
        f"filtered/untargeted audit rows: {len(rtk_filtered_out)}."
    )

    keep_cols = (
        EXPECTED_CLASSIFIER_COLS
        + OPTIONAL_CONTEXT_COLS
        + [
            ACTUAL_DIRECTION_COL_NAME,
            "_composite_point_key",
            "_is_in_anomaly",
            "classifier_file",
            "classifier_file_path",
            "classifier_id",
            "evaluation_method",
            "evaluation_method_analysis",
            "rtk_lookup_source",
            "rtk_pair_level_fallback",
            "rtk_pair_level_valid_compare_to_rows",
            "rtk_pair_level_difference_rows",
            "rtk_pair_level_same_rows",
            "rtk_raw_compare_to_values",
            "rtk_normalized_compare_to_values",
            "rtk_pairwise_evidence_rows",
            "rtk_pairwise_evaluable_rows",
            "rtk_pairwise_distinct_valid_directions",
            "rtk_pairwise_inverted_rows",
            "rtk_pairwise_role_counts",
            "rtk_pairwise_source_keys",
            "rtk_distinct_stads_directions",
            "rtk_stads_direction_values",
            "rtk_joined_anomaly_ids",
            "rtk_analysis_key_was_made_anomaly_specific",
        ]
    )
    keep_cols = [col for col in dict.fromkeys(keep_cols) if col in rtk_in.columns]

    if not pairwise_evidence.empty:
        pairwise_evidence = pairwise_evidence.copy()
        pairwise_evidence["rtk_lookup_source"] = lookup_source
    else:
        pairwise_evidence = key_overlap_summary.copy()
        pairwise_evidence["rtk_lookup_source"] = lookup_source
        pairwise_evidence["rtk_pairwise_status"] = "no_pairwise_evidence_created"

    return rtk_in[keep_cols].copy(), rtk_filtered_out.copy(), pairwise_evidence.copy()

def clean_and_derive_rating_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize reviewer ratings and derive analysis columns."""
    df = df.copy()

    cols_to_clean = [
        "Crop",
        "Anomaly Indicated",
        "Anomaly Classification Confidence",
        "Direction Indicated",
        "Direction Confidence",
        ACTUAL_DIRECTION_COL_NAME,
    ]

    for col in cols_to_clean:
        if col in df.columns:
            df[col] = clean_text_col(df[col])

    if "Crop" in df.columns:
        df["Crop"] = df["Crop"].replace(
            {
                "Corn ": "Corn",
                "Soybeans": "Soybean",
                "Soybean ": "Soybean",
            }
        )

    df["anomaly_indicated_normalized"] = df["Anomaly Indicated"].apply(
        normalize_anomaly_value
    )
    df["anomaly_confidence_normalized"] = df[
        "Anomaly Classification Confidence"
    ].apply(normalize_confidence_value)
    df["direction_indicated_normalized"] = df["Direction Indicated"].apply(
        normalize_direction_value
    )
    df["direction_confidence_normalized"] = df["Direction Confidence"].apply(
        normalize_confidence_value
    )
    df["stads_direction_normalized"] = df[ACTUAL_DIRECTION_COL_NAME].apply(
        normalize_direction_value
    )

    df["anomaly_support_score"] = df["anomaly_indicated_normalized"].map(
        ANOMALY_SUPPORT_SCORE
    )
    df["anomaly_confidence_weight"] = df["anomaly_confidence_normalized"].map(
        CONFIDENCE_WEIGHT
    )

    # Rating-level binary confirmation outcome for performance modeling.
    df["confirmation_binary_moderate"] = df["anomaly_indicated_normalized"].map(
        {
            "Confirmed": 1,
            "Possibly Confirmed": 1,
            "Not Confirmed": 0,
        }
    )

    # Rating-level direction match status.
    valid_direction_mask = (
        df["direction_indicated_normalized"].isin(["Positive", "Negative"])
        & df["stads_direction_normalized"].isin(["Positive", "Negative"])
    )
    df["direction_match_status"] = "Indeterminate"
    df.loc[
        valid_direction_mask
        & (df["direction_indicated_normalized"] == df["stads_direction_normalized"]),
        "direction_match_status",
    ] = "Match"
    df.loc[
        valid_direction_mask
        & (df["direction_indicated_normalized"] != df["stads_direction_normalized"]),
        "direction_match_status",
    ] = "Mismatch"
    df["direction_match_binary"] = df["direction_match_status"].map(
        {"Match": True, "Mismatch": False}
    )

    # Evaluation-method detection. If absent, default to Blind paired because the
    # current file list uses BlindAnomalyClassification_* names.
    method_col = detect_evaluation_method_column(df)
    if method_col:
        df["evaluation_method_analysis"] = clean_text_col(df[method_col])
    else:
        df["evaluation_method_analysis"] = "Blind paired"

    return df


def detect_evaluation_method_column(df: pd.DataFrame) -> Optional[str]:
    """Return the first available evaluation-method column, if any."""
    for col in EVALUATION_METHOD_COL_CANDIDATES:
        if col in df.columns:
            return col
    return None


# =============================================================================
# INTER-RATER RELIABILITY
# =============================================================================


def make_rater_matrix(
    df: pd.DataFrame,
    rating_col: str,
    item_col: str = "_composite_point_key",
    rater_col: str = "classifier_id",
    valid_values: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Create item x rater matrix for reliability analysis."""
    tmp = df[[item_col, rater_col, rating_col]].copy()
    tmp = tmp.dropna(subset=[item_col, rater_col, rating_col])

    if valid_values is not None:
        tmp = tmp[tmp[rating_col].isin(valid_values)].copy()

    # If a reviewer has multiple ratings for the same point, use the unique mode.
    # Ties are set to missing because they cannot be represented as one rating.
    tmp = (
        tmp.groupby([item_col, rater_col], dropna=False)[rating_col]
        .agg(mode_or_tie)
        .reset_index()
    )
    tmp.loc[tmp[rating_col] == "Tie", rating_col] = pd.NA

    matrix = tmp.pivot(index=item_col, columns=rater_col, values=rating_col)
    return matrix


def build_coincidence_matrix(
    ratings_matrix: pd.DataFrame,
    categories: list[str],
) -> pd.DataFrame:
    """Build Krippendorff coincidence matrix."""
    coincidence = pd.DataFrame(0.0, index=categories, columns=categories)

    for _, row in ratings_matrix.iterrows():
        values = row.dropna().tolist()
        values = [v for v in values if v in categories]
        m = len(values)
        if m <= 1:
            continue

        counts = pd.Series(values).value_counts()
        for c in categories:
            for k in categories:
                n_c = counts.get(c, 0)
                n_k = counts.get(k, 0)
                if c == k:
                    coincidence.loc[c, k] += n_c * (n_c - 1) / (m - 1)
                else:
                    coincidence.loc[c, k] += n_c * n_k / (m - 1)

    return coincidence


def krippendorff_alpha(
    ratings_matrix: pd.DataFrame,
    categories: Optional[list[str]] = None,
    distance: str = "nominal",
    numeric_scores: Optional[dict[str, float]] = None,
) -> float:
    """Compute Krippendorff's alpha from an item x rater matrix.

    distance='nominal' uses 0/1 disagreement.
    distance='interval' uses squared distance among numeric_scores.
    """
    if categories is None:
        categories = sorted(pd.unique(ratings_matrix.values.ravel()))
        categories = [c for c in categories if pd.notna(c)]

    categories = [c for c in categories if c in set(pd.Series(ratings_matrix.values.ravel()).dropna())]
    if len(categories) <= 1:
        return np.nan

    coincidence = build_coincidence_matrix(ratings_matrix, categories)
    N = coincidence.to_numpy().sum()
    if N <= 1:
        return np.nan

    if distance == "nominal":
        # Build the distance matrix as a writable NumPy array first.
        # On some pandas/NumPy versions, `delta.values` from a DataFrame
        # can be read-only, which causes np.fill_diagonal(...) to fail.
        delta_array = np.ones((len(categories), len(categories)), dtype=float)
        np.fill_diagonal(delta_array, 0.0)
        delta = pd.DataFrame(delta_array, index=categories, columns=categories)
    elif distance == "interval":
        if numeric_scores is None:
            raise ValueError("numeric_scores is required for interval Krippendorff alpha")
        delta = pd.DataFrame(index=categories, columns=categories, dtype=float)
        for c in categories:
            for k in categories:
                delta.loc[c, k] = (numeric_scores[c] - numeric_scores[k]) ** 2
    else:
        raise ValueError("distance must be 'nominal' or 'interval'")

    observed_disagreement = (coincidence * delta).to_numpy().sum() / N

    marginals = coincidence.sum(axis=1)
    expected_disagreement = 0.0
    for c in categories:
        for k in categories:
            expected_disagreement += marginals[c] * marginals[k] * delta.loc[c, k]
    expected_disagreement = expected_disagreement / (N * (N - 1))

    if expected_disagreement == 0:
        return np.nan
    return float(1 - observed_disagreement / expected_disagreement)


def fleiss_kappa_from_matrix(
    ratings_matrix: pd.DataFrame,
    categories: Optional[list[str]] = None,
) -> float:
    """Compute Fleiss' kappa on complete rows with the same number of ratings."""
    row_counts = ratings_matrix.notna().sum(axis=1)
    if row_counts.empty or row_counts.max() < 2:
        return np.nan

    target_n = int(row_counts.mode().iloc[0])
    complete = ratings_matrix[row_counts == target_n].copy()
    if complete.empty:
        return np.nan

    if categories is None:
        categories = sorted(pd.unique(complete.values.ravel()))
        categories = [c for c in categories if pd.notna(c)]
    if len(categories) <= 1:
        return np.nan

    count_matrix = []
    for _, row in complete.iterrows():
        counts = row.value_counts()
        count_matrix.append([counts.get(cat, 0) for cat in categories])

    count_matrix = np.asarray(count_matrix, dtype=float)
    N, k = count_matrix.shape
    n = count_matrix.sum(axis=1)
    if N == 0 or np.any(n < 2) or len(set(n)) != 1:
        return np.nan

    n_raters = n[0]
    P_i = ((count_matrix**2).sum(axis=1) - n_raters) / (n_raters * (n_raters - 1))
    P_bar = P_i.mean()
    p_j = count_matrix.sum(axis=0) / (N * n_raters)
    P_e = (p_j**2).sum()

    if P_e == 1:
        return np.nan
    return float((P_bar - P_e) / (1 - P_e))


def gwet_ac1_from_matrix(
    ratings_matrix: pd.DataFrame,
    categories: Optional[list[str]] = None,
) -> float:
    """Compute nominal Gwet's AC1 for multiple raters."""
    if categories is None:
        categories = sorted(pd.unique(ratings_matrix.values.ravel()))
        categories = [c for c in categories if pd.notna(c)]
    categories = [c for c in categories if pd.notna(c)]
    q = len(categories)
    if q <= 1:
        return np.nan

    observed_agreements = []
    all_values = []

    for _, row in ratings_matrix.iterrows():
        values = [v for v in row.dropna().tolist() if v in categories]
        m = len(values)
        if m <= 1:
            continue
        counts = pd.Series(values).value_counts()
        agreement = sum(count * (count - 1) for count in counts) / (m * (m - 1))
        observed_agreements.append(agreement)
        all_values.extend(values)

    if not observed_agreements or not all_values:
        return np.nan

    P_a = float(np.mean(observed_agreements))
    p = pd.Series(all_values).value_counts(normalize=True).reindex(categories, fill_value=0)
    P_e = float((p * (1 - p)).sum() / (q - 1))

    if P_e == 1:
        return np.nan
    return float((P_a - P_e) / (1 - P_e))


def cohen_kappa(y1: pd.Series, y2: pd.Series, categories: Optional[list[str]] = None) -> float:
    """Compute unweighted Cohen's kappa for two raters."""
    paired = pd.DataFrame({"y1": y1, "y2": y2}).dropna()
    if paired.empty:
        return np.nan

    if categories is None:
        categories = sorted(set(paired["y1"]).union(set(paired["y2"])))
    if len(categories) <= 1:
        return np.nan

    obs = (paired["y1"] == paired["y2"]).mean()
    p1 = paired["y1"].value_counts(normalize=True).reindex(categories, fill_value=0)
    p2 = paired["y2"].value_counts(normalize=True).reindex(categories, fill_value=0)
    exp = float((p1 * p2).sum())

    if exp == 1:
        return np.nan
    return float((obs - exp) / (1 - exp))


def pairwise_kappas(ratings_matrix: pd.DataFrame, categories: Optional[list[str]] = None) -> pd.DataFrame:
    """Compute pairwise Cohen's kappa among all reviewer columns."""
    rows = []
    for r1, r2 in itertools.combinations(ratings_matrix.columns, 2):
        paired = ratings_matrix[[r1, r2]].dropna()
        rows.append(
            {
                "rater_1": r1,
                "rater_2": r2,
                "n_common_items": len(paired),
                "raw_agreement": (paired[r1] == paired[r2]).mean() if len(paired) else np.nan,
                "cohen_kappa": cohen_kappa(paired[r1], paired[r2], categories),
            }
        )
    return pd.DataFrame(rows)


def run_inter_rater_reliability(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Run all reliability analyses and return output tables."""
    outputs: dict[str, pd.DataFrame] = {}

    reliability_specs = [
        {
            "analysis": "anomaly_all_categories",
            "rating_col": "anomaly_indicated_normalized",
            "valid_values": ANOMALY_INDICATED_ORDER,
            "nominal_categories": ANOMALY_INDICATED_ORDER,
            "interval_scores": None,
        },
        {
            "analysis": "anomaly_evaluable_only",
            "rating_col": "anomaly_indicated_normalized",
            "valid_values": ["Confirmed", "Possibly Confirmed", "Not Confirmed"],
            "nominal_categories": ["Confirmed", "Possibly Confirmed", "Not Confirmed"],
            "interval_scores": ANOMALY_SUPPORT_SCORE,
        },
        {
            "analysis": "direction_evaluable_only",
            "rating_col": "direction_indicated_normalized",
            "valid_values": ["Positive", "Negative", "Same"],
            "nominal_categories": ["Positive", "Negative", "Same"],
            "interval_scores": None,
        },
    ]

    summary_rows = []

    for spec in reliability_specs:
        matrix = make_rater_matrix(
            df,
            rating_col=spec["rating_col"],
            valid_values=spec["valid_values"],
        )
        categories = spec["nominal_categories"]

        summary_rows.append(
            {
                "analysis": spec["analysis"],
                "n_items_with_any_rating": int(matrix.notna().any(axis=1).sum()),
                "n_items_with_2plus_ratings": int((matrix.notna().sum(axis=1) >= 2).sum()),
                "n_raters": int(matrix.shape[1]),
                "krippendorff_alpha_nominal": krippendorff_alpha(
                    matrix, categories=categories, distance="nominal"
                ),
                "krippendorff_alpha_interval": (
                    krippendorff_alpha(
                        matrix,
                        categories=list(spec["interval_scores"].keys()),
                        distance="interval",
                        numeric_scores=spec["interval_scores"],
                    )
                    if spec["interval_scores"] is not None
                    else np.nan
                ),
                "fleiss_kappa": fleiss_kappa_from_matrix(matrix, categories=categories),
                "gwet_ac1": gwet_ac1_from_matrix(matrix, categories=categories),
            }
        )

        outputs[f"rater_matrix_{spec['analysis']}"] = matrix.reset_index()
        outputs[f"pairwise_kappa_{spec['analysis']}"] = pairwise_kappas(
            matrix, categories=categories
        )

    outputs["reliability_summary"] = pd.DataFrame(summary_rows)
    return outputs


# =============================================================================
# CONSENSUS AND PERFORMANCE SUMMARIES
# =============================================================================


def build_point_level_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse reviewer ratings to one consensus row per reviewed point."""
    rows = []

    metadata_cols = [
        "Crop",
        "fn",
        "farmnumber",
        "fieldname",
        "sent_lon",
        "sent_lat",
        "stads_direction_normalized",
        ACTUAL_DIRECTION_COL_NAME,
        "evaluation_method_analysis",
        # RTK joined-anomaly and issue metadata retained at point/anomaly level
        "sample_pai",
        "a_or_b",
        "compare_to",
        "bad_things",
        "bad_thin_1",
        "crop_obser",
        "rtk_lookup_source",
        "rtk_joined_to_anomaly",
        "rtk_joined_anomaly_id",
        "rtk_joined_anomaly_ids",
        "rtk_joined_anomaly_point_id",
        "rtk_distinct_stads_directions",
        "rtk_stads_direction_values",
        "CropName",
        "imagerydate",
        "layer",
        "path",
    ]
    metadata_cols = [col for col in metadata_cols if col in df.columns]

    numeric_optional_cols = [
        col
        for col in OPTIONAL_CONTEXT_COLS
        if col in df.columns and pd.api.types.is_numeric_dtype(df[col])
    ]

    for point_key, group in df.groupby("_composite_point_key", dropna=False):
        n_ratings_total = len(group)
        anomaly_valid = group[
            group["anomaly_indicated_normalized"].isin(
                ["Confirmed", "Possibly Confirmed", "Not Confirmed"]
            )
        ].copy()
        n_anomaly_valid = len(anomaly_valid)

        n_confirmed = int((anomaly_valid["anomaly_indicated_normalized"] == "Confirmed").sum())
        n_possible = int((anomaly_valid["anomaly_indicated_normalized"] == "Possibly Confirmed").sum())
        n_not = int((anomaly_valid["anomaly_indicated_normalized"] == "Not Confirmed").sum())
        n_supporting = n_confirmed + n_possible

        confidence_weights = anomaly_valid["anomaly_confidence_weight"]
        support_scores = anomaly_valid["anomaly_support_score"]
        valid_weight_mask = confidence_weights.notna() & support_scores.notna()
        confidence_weighted_support = (
            np.average(
                support_scores[valid_weight_mask].astype(float),
                weights=confidence_weights[valid_weight_mask].astype(float),
            )
            if valid_weight_mask.any()
            else np.nan
        )

        mean_support_score = (
            float(anomaly_valid["anomaly_support_score"].mean())
            if n_anomaly_valid > 0
            else np.nan
        )

        evaluation_method = (
            first_non_missing(group["evaluation_method_analysis"])
            if "evaluation_method_analysis" in group.columns
            else BLIND_METHOD_LABEL
        )
        min_raters_for_point = 1 if evaluation_method == RTK_METHOD_LABEL else MIN_RATERS_FOR_CONSENSUS

        direction_valid = group[
            group["direction_indicated_normalized"].isin(["Positive", "Negative", "Same"])
        ].copy()
        n_direction_valid = len(direction_valid)
        consensus_direction = mode_or_tie(direction_valid["direction_indicated_normalized"])

        stads_direction = first_non_missing(group["stads_direction_normalized"])
        direction_consensus_match = pd.NA
        if (
            pd.notna(consensus_direction)
            and pd.notna(stads_direction)
            and consensus_direction in ["Positive", "Negative"]
            and stads_direction in ["Positive", "Negative"]
        ):
            direction_consensus_match = consensus_direction == stads_direction

        row = {
            "_composite_point_key": point_key,
            "n_ratings_total": n_ratings_total,
            "n_anomaly_valid": n_anomaly_valid,
            "n_confirmed": n_confirmed,
            "n_possibly_confirmed": n_possible,
            "n_not_confirmed": n_not,
            "n_supporting_confirmed_or_possible": n_supporting,
            "evaluation_method_analysis": evaluation_method,
            "min_raters_required_for_rule": min_raters_for_point,
            "mean_support_score": mean_support_score,
            "confidence_weighted_support": confidence_weighted_support,
            "strict_confirmed_3of4_or_75pct": (
                n_anomaly_valid >= min_raters_for_point
                and n_confirmed >= threshold_count(n_anomaly_valid, 0.75)
            ),
            "moderate_confirmed_3of4_or_75pct_confirmed_possible": (
                n_anomaly_valid >= min_raters_for_point
                and n_supporting >= threshold_count(n_anomaly_valid, 0.75)
            ),
            "liberal_confirmed_2of4_or_50pct_confirmed_possible": (
                n_anomaly_valid >= min_raters_for_point
                and n_supporting >= threshold_count(n_anomaly_valid, 0.50)
            ),
            "confidence_weighted_confirmed": (
                pd.notna(confidence_weighted_support)
                and confidence_weighted_support >= CONFIDENCE_WEIGHTED_SUPPORT_THRESHOLD
            ),
            "consensus_anomaly_label_mode": mode_or_tie(
                anomaly_valid["anomaly_indicated_normalized"]
            ),
            "n_direction_valid": n_direction_valid,
            "consensus_direction_mode": consensus_direction,
            "stads_direction_normalized": stads_direction,
            "direction_consensus_match": direction_consensus_match,
        }

        for col in metadata_cols:
            if col not in row:
                row[col] = first_non_missing(group[col])

        for col in numeric_optional_cols:
            row[col] = pd.to_numeric(group[col], errors="coerce").mean()

        rows.append(row)

    consensus = pd.DataFrame(rows)
    return consensus


def run_consensus_support_analysis(consensus: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Summarize STADS support under multiple consensus rules."""
    outputs: dict[str, pd.DataFrame] = {}

    rule_cols = [
        "strict_confirmed_3of4_or_75pct",
        "moderate_confirmed_3of4_or_75pct_confirmed_possible",
        "liberal_confirmed_2of4_or_50pct_confirmed_possible",
        "confidence_weighted_confirmed",
    ]

    rows = []
    if "min_raters_required_for_rule" in consensus.columns:
        evaluable = consensus[
            consensus["n_anomaly_valid"] >= consensus["min_raters_required_for_rule"]
        ].copy()
    else:
        evaluable = consensus[consensus["n_anomaly_valid"] >= MIN_RATERS_FOR_CONSENSUS].copy()
    for rule_col in rule_cols:
        rate_table = summarize_binary_rate(evaluable, rule_col, label="confirmation_rate")
        row = rate_table.iloc[0].to_dict()
        row["consensus_rule"] = rule_col
        rows.append(row)

    outputs["consensus_rule_sensitivity"] = pd.DataFrame(rows)[
        [
            "consensus_rule",
            "n",
            "successes",
            "confirmation_rate",
            "confirmation_rate_percent",
            "ci_lower",
            "ci_upper",
            "ci_lower_percent",
            "ci_upper_percent",
        ]
    ]

    if "Crop" in consensus.columns:
        crop_rows = []
        for rule_col in rule_cols:
            tmp = summarize_binary_rate(
                evaluable,
                rule_col,
                group_cols=["Crop"],
                label="confirmation_rate",
            )
            tmp["consensus_rule"] = rule_col
            crop_rows.append(tmp)
        outputs["consensus_by_crop"] = pd.concat(crop_rows, ignore_index=True)

    if "stads_direction_normalized" in consensus.columns:
        direction_rows = []
        for rule_col in rule_cols:
            tmp = summarize_binary_rate(
                evaluable,
                rule_col,
                group_cols=["stads_direction_normalized"],
                label="confirmation_rate",
            )
            tmp["consensus_rule"] = rule_col
            direction_rows.append(tmp)
        outputs["consensus_by_stads_direction"] = pd.concat(direction_rows, ignore_index=True)

    # Distribution of continuous support scores.
    outputs["support_score_summary"] = (
        evaluable[["mean_support_score", "confidence_weighted_support"]]
        .describe()
        .reset_index()
        .rename(columns={"index": "statistic"})
    )

    return outputs


def run_direction_agreement_analysis(
    df: pd.DataFrame,
    consensus: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Summarize reviewer-level and consensus-level direction agreement."""
    outputs: dict[str, pd.DataFrame] = {}

    rating_evaluable = df[df["direction_match_binary"].notna()].copy()
    outputs["direction_rating_overall"] = summarize_binary_rate(
        rating_evaluable,
        "direction_match_binary",
        label="direction_agreement_rate",
    )
    outputs["direction_rating_by_crop"] = summarize_binary_rate(
        rating_evaluable,
        "direction_match_binary",
        group_cols=["Crop"] if "Crop" in rating_evaluable.columns else None,
        label="direction_agreement_rate",
    )
    outputs["direction_rating_by_stads_dir"] = summarize_binary_rate(
        rating_evaluable,
        "direction_match_binary",
        group_cols=["stads_direction_normalized"],
        label="direction_agreement_rate",
    )
    if "evaluation_method_analysis" in rating_evaluable.columns:
        outputs["direction_rating_by_method"] = summarize_binary_rate(
            rating_evaluable,
            "direction_match_binary",
            group_cols=["evaluation_method_analysis"],
            label="direction_agreement_rate",
        )
        if "Crop" in rating_evaluable.columns:
            outputs["direction_rating_by_method_crop"] = summarize_binary_rate(
                rating_evaluable,
                "direction_match_binary",
                group_cols=["evaluation_method_analysis", "Crop"],
                label="direction_agreement_rate",
            )

    consensus_evaluable = consensus[consensus["direction_consensus_match"].notna()].copy()
    if not consensus_evaluable.empty:
        consensus_evaluable["direction_consensus_match"] = consensus_evaluable[
            "direction_consensus_match"
        ].astype(bool)
        outputs["direction_consensus_overall"] = summarize_binary_rate(
            consensus_evaluable,
            "direction_consensus_match",
            label="direction_agreement_rate",
        )
        if "Crop" in consensus_evaluable.columns:
            outputs["direction_consensus_by_crop"] = summarize_binary_rate(
                consensus_evaluable,
                "direction_consensus_match",
                group_cols=["Crop"],
                label="direction_agreement_rate",
            )
        outputs["direction_consensus_by_stads_dir"] = summarize_binary_rate(
            consensus_evaluable,
            "direction_consensus_match",
            group_cols=["stads_direction_normalized"],
            label="direction_agreement_rate",
        )
        if "evaluation_method_analysis" in consensus_evaluable.columns:
            outputs["direction_consensus_by_method"] = summarize_binary_rate(
                consensus_evaluable,
                "direction_consensus_match",
                group_cols=["evaluation_method_analysis"],
                label="direction_agreement_rate",
            )
            if "Crop" in consensus_evaluable.columns:
                outputs["direction_consensus_by_method_crop"] = summarize_binary_rate(
                    consensus_evaluable,
                    "direction_consensus_match",
                    group_cols=["evaluation_method_analysis", "Crop"],
                    label="direction_agreement_rate",
                )

        # One-sided binomial test against 0.5 for overall consensus direction agreement.
        k = int(consensus_evaluable["direction_consensus_match"].sum())
        n = int(len(consensus_evaluable))
        outputs["direction_consensus_binomial_test"] = pd.DataFrame(
            [
                {
                    "successes": k,
                    "n": n,
                    "null_probability": 0.5,
                    "one_sided_p_value_greater_than_chance": exact_binomial_p_value_greater_equal(k, n, 0.5),
                }
            ]
        )
    else:
        outputs["direction_consensus_overall"] = pd.DataFrame(
            [{"note": "No evaluable consensus direction rows."}]
        )

    return outputs


def rename_consensus_outputs(
    outputs: dict[str, pd.DataFrame],
    prefix: str,
) -> dict[str, pd.DataFrame]:
    """Rename generic consensus-support outputs using a method prefix.

    This keeps the workbook/table names parallel across evidence streams. For
    example, the generic consensus_rule_sensitivity table becomes
    blind_confirmation_rates or rtk_confirmation_rates depending on prefix.
    """
    name_map = {
        "consensus_rule_sensitivity": f"{prefix}_confirmation_rates",
        "consensus_by_crop": f"{prefix}_confirmation_by_crop",
        "consensus_by_stads_direction": f"{prefix}_confirmation_by_stads_direction",
        "support_score_summary": f"{prefix}_support_score_summary",
    }
    return {name_map.get(name, f"{prefix}_{name}"): table for name, table in outputs.items()}


def rename_direction_outputs(
    outputs: dict[str, pd.DataFrame],
    prefix: str,
) -> dict[str, pd.DataFrame]:
    """Rename generic direction-agreement outputs using a method prefix.

    The consensus-level direction-agreement table receives the simplest name,
    such as blind_direction_agreement or rtk_direction_agreement, because this is
    the main direction-performance summary for each evidence stream.
    """
    name_map = {
        "direction_rating_overall": f"{prefix}_direction_rating_agreement",
        "direction_rating_by_crop": f"{prefix}_direction_rating_by_crop",
        "direction_rating_by_stads_dir": f"{prefix}_direction_rating_by_stads_direction",
        "direction_rating_by_method": f"{prefix}_direction_rating_by_method",
        "direction_rating_by_method_crop": f"{prefix}_direction_rating_by_method_crop",
        "direction_consensus_overall": f"{prefix}_direction_agreement",
        "direction_consensus_by_crop": f"{prefix}_direction_agreement_by_crop",
        "direction_consensus_by_stads_dir": f"{prefix}_direction_agreement_by_stads_direction",
        "direction_consensus_by_method": f"{prefix}_direction_agreement_by_method",
        "direction_consensus_by_method_crop": f"{prefix}_direction_agreement_by_method_crop",
        "direction_consensus_binomial_test": f"{prefix}_direction_binomial_test",
    }
    return {name_map.get(name, f"{prefix}_{name}"): table for name, table in outputs.items()}


def make_parallel_method_summary_index(
    blind_consensus_outputs: dict[str, pd.DataFrame],
    rtk_consensus_outputs: dict[str, pd.DataFrame],
    blind_direction_outputs: dict[str, pd.DataFrame],
    rtk_direction_outputs: dict[str, pd.DataFrame],
    method_outputs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Create an index explaining the main parallel output tables."""
    expected = [
        ("blind_confirmation_rates", "Blind paired confirmation-rate sensitivity across consensus rules."),
        ("rtk_confirmation_rates", "RTK-guided confirmation-rate sensitivity across the same rule columns; each RTK point is a single-observer evidence row."),
        ("method_confirmation_rates", "Direct method comparison using the moderate confirmation rule."),
        ("blind_direction_agreement", "Blind paired consensus-level direction agreement with STADS direction."),
        ("rtk_direction_agreement", "RTK-guided point-level direction agreement with STADS direction."),
        ("blind_confirmation_by_crop", "Blind paired confirmation rates by crop."),
        ("rtk_confirmation_by_crop", "RTK-guided confirmation rates by crop."),
        ("blind_direction_agreement_by_crop", "Blind paired direction agreement by crop."),
        ("rtk_direction_agreement_by_crop", "RTK-guided direction agreement by crop."),
    ]
    available = set(blind_consensus_outputs) | set(rtk_consensus_outputs) | set(blind_direction_outputs) | set(rtk_direction_outputs) | set(method_outputs)
    return pd.DataFrame(
        [
            {
                "table_name": name,
                "available_in_workbook": name in available,
                "description": description,
            }
            for name, description in expected
        ]
    )


# =============================================================================
# MIXED-EFFECTS / CLUSTERED LOGISTIC MODEL
# =============================================================================


def choose_model_predictors(model_df: pd.DataFrame) -> list[str]:
    """Choose available fixed-effect predictors for the confirmation model."""
    predictors = []

    categorical_candidates = ["Crop", "stads_direction_normalized"]
    for col in categorical_candidates:
        if col in model_df.columns and model_df[col].nunique(dropna=True) >= 2:
            predictors.append(f"C({col})")

    numeric_candidates = [
        "anomaly_area_m2",
        "anomaly_area_ha",
        "area_m2",
        "area_ha",
        "stads_magnitude",
        "anomaly_magnitude",
        "magnitude",
        "ndvi_difference",
        "ndvi_diff",
        "mean_deviation",
        "max_deviation",
        "persistence",
        "days_detected",
    ]
    for col in numeric_candidates:
        if col in model_df.columns:
            numeric = pd.to_numeric(model_df[col], errors="coerce")
            if numeric.notna().sum() >= 10 and numeric.nunique(dropna=True) >= 2:
                model_df[col] = numeric
                predictors.append(col)

    if not predictors:
        predictors = ["1"]

    return predictors


def run_mixed_effects_model(df: pd.DataFrame) -> dict[str, object]:
    """Fit optional mixed / clustered logistic model and return model outputs."""
    outputs: dict[str, object] = {}

    model_df = df[df["confirmation_binary_moderate"].notna()].copy()
    model_df["confirmation_binary_moderate"] = model_df[
        "confirmation_binary_moderate"
    ].astype(int)

    if len(model_df) < 20 or model_df["confirmation_binary_moderate"].nunique() < 2:
        outputs["model_note"] = (
            "Model not fit: fewer than 20 evaluable reviewer ratings or no variation "
            "in confirmation_binary_moderate."
        )
        return outputs

    predictors = choose_model_predictors(model_df)
    formula = "confirmation_binary_moderate ~ " + " + ".join(predictors)
    outputs["model_formula"] = formula
    outputs["n_model_rows"] = len(model_df)

    # First try a Bayesian generalized linear mixed model with reviewer and field
    # random effects. This is closer to the desired mixed-effects model.
    try:
        from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

        vc_formulas = {}
        if "classifier_id" in model_df.columns and model_df["classifier_id"].nunique() >= 2:
            vc_formulas["reviewer"] = "0 + C(classifier_id)"
        if "fn" in model_df.columns and model_df["fn"].nunique(dropna=True) >= 2:
            vc_formulas["field"] = "0 + C(fn)"

        if vc_formulas:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = BinomialBayesMixedGLM.from_formula(formula, vc_formulas, model_df)
                fit = model.fit_vb()
            outputs["mixed_model_type"] = "BinomialBayesMixedGLM.fit_vb"
            outputs["mixed_model_summary_text"] = str(fit.summary())
            return outputs
    except Exception as exc:
        outputs["mixed_model_error"] = repr(exc)

    # Fallback: logistic GLM with reviewer fixed effect and clustered standard
    # errors by point key. This is not a true mixed model, but it is useful and
    # robust for publication diagnostics.
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        fallback_formula = formula
        if "classifier_id" in model_df.columns and model_df["classifier_id"].nunique() >= 2:
            fallback_formula += " + C(classifier_id)"

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            glm_fit = smf.glm(
                fallback_formula,
                data=model_df,
                family=sm.families.Binomial(),
            ).fit(
                cov_type="cluster",
                cov_kwds={"groups": model_df["_composite_point_key"]},
            )

        outputs["mixed_model_type"] = "GLM Binomial with reviewer fixed effects and clustered SE by point"
        outputs["mixed_model_summary_text"] = str(glm_fit.summary())
    except Exception as exc:
        outputs["glm_fallback_error"] = repr(exc)
        outputs["model_note"] = (
            "Model could not be fit. Install statsmodels or simplify the predictor set."
        )

    return outputs


# =============================================================================
# BLIND VS RTK-GUIDED COMPARISON
# =============================================================================


def run_method_comparison(consensus: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Compare confirmation rates between blind and RTK-guided evidence streams."""
    outputs: dict[str, pd.DataFrame] = {}

    if "evaluation_method_analysis" not in consensus.columns:
        outputs["method_comparison_note"] = pd.DataFrame(
            [
                {
                    "note": (
                        "No evaluation-method column found. Add a column such as "
                        "evaluation_method or scouting_method with values like "
                        "'Blind paired' and 'RTK guided' to run this comparison."
                    )
                }
            ]
        )
        return outputs

    n_methods = consensus["evaluation_method_analysis"].nunique(dropna=True)
    if n_methods < 2:
        outputs["method_comparison_note"] = pd.DataFrame(
            [
                {
                    "note": (
                        "Only one evaluation method was detected: "
                        f"{consensus['evaluation_method_analysis'].dropna().unique().tolist()}. "
                        "Blind vs RTK-guided comparison requires both evidence streams "
                        "or a method/source column in the input data."
                    )
                }
            ]
        )
        return outputs

    if "min_raters_required_for_rule" in consensus.columns:
        evaluable = consensus[
            consensus["n_anomaly_valid"] >= consensus["min_raters_required_for_rule"]
        ].copy()
    else:
        evaluable = consensus[consensus["n_anomaly_valid"] >= MIN_RATERS_FOR_CONSENSUS].copy()
    method_rows = []
    # Use the moderate unweighted rule for method comparison. RTK-guided
    # assessments do not have independent reviewer confidence scores, so
    # confidence-weighted consensus is intentionally not compared across methods.
    for rule_col in [
        "moderate_confirmed_3of4_or_75pct_confirmed_possible",
    ]:
        tmp = summarize_binary_rate(
            evaluable,
            rule_col,
            group_cols=["evaluation_method_analysis"],
            label="confirmation_rate",
        )
        tmp["consensus_rule"] = rule_col
        method_rows.append(tmp)

    outputs["method_confirmation_rates"] = pd.concat(method_rows, ignore_index=True)

    # Simple model-based comparison if statsmodels is available.
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        model_df = evaluable[
            [
                "moderate_confirmed_3of4_or_75pct_confirmed_possible",
                "evaluation_method_analysis",
                "Crop" if "Crop" in evaluable.columns else "evaluation_method_analysis",
            ]
        ].copy()
        model_df = model_df.rename(
            columns={"moderate_confirmed_3of4_or_75pct_confirmed_possible": "confirmed"}
        )
        model_df["confirmed"] = model_df["confirmed"].astype(int)

        formula = "confirmed ~ C(evaluation_method_analysis)"
        if "Crop" in evaluable.columns and evaluable["Crop"].nunique(dropna=True) >= 2:
            model_df["Crop"] = evaluable["Crop"]
            formula += " + C(Crop)"

        if model_df["confirmed"].nunique() >= 2:
            fit = smf.glm(formula, data=model_df, family=sm.families.Binomial()).fit()
            outputs["method_logistic_model_summary"] = pd.DataFrame(
                {"model_summary": str(fit.summary()).splitlines()}
            )
    except Exception as exc:
        outputs["method_model_note"] = pd.DataFrame(
            [{"note": f"Method-comparison model not fit: {repr(exc)}"}]
        )

    return outputs


# =============================================================================
# CONFIDENCE ANALYSIS
# =============================================================================


def run_confidence_analysis(df: pd.DataFrame, consensus: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Analyze whether reviewer confidence corresponds to consensus/agreement."""
    outputs: dict[str, pd.DataFrame] = {}

    consensus_cols = [
        "_composite_point_key",
        "consensus_anomaly_label_mode",
        "moderate_confirmed_3of4_or_75pct_confirmed_possible",
        "confidence_weighted_support",
    ]
    merged = df.merge(
        consensus[consensus_cols],
        on="_composite_point_key",
        how="left",
        suffixes=("", "_consensus"),
    )

    merged["matches_consensus_anomaly_label"] = pd.NA
    mask = (
        merged["anomaly_indicated_normalized"].isin(
            ["Confirmed", "Possibly Confirmed", "Not Confirmed"]
        )
        & merged["consensus_anomaly_label_mode"].isin(
            ["Confirmed", "Possibly Confirmed", "Not Confirmed"]
        )
    )
    merged.loc[mask, "matches_consensus_anomaly_label"] = (
        merged.loc[mask, "anomaly_indicated_normalized"]
        == merged.loc[mask, "consensus_anomaly_label_mode"]
    )

    # Confidence vs match to consensus anomaly label.
    confidence_match = merged[merged["matches_consensus_anomaly_label"].notna()].copy()
    confidence_match["matches_consensus_anomaly_label"] = confidence_match[
        "matches_consensus_anomaly_label"
    ].astype(bool)
    if not confidence_match.empty:
        outputs["anomaly_confidence_vs_consensus"] = summarize_binary_rate(
            confidence_match,
            "matches_consensus_anomaly_label",
            group_cols=["anomaly_confidence_normalized"],
            label="consensus_agreement_rate",
        )

    # Confidence vs rating-level confirmed/possible rate.
    anomaly_confidence_summary = (
        merged[merged["confirmation_binary_moderate"].notna()]
        .groupby("anomaly_confidence_normalized", dropna=False)
        .agg(
            n=("confirmation_binary_moderate", "size"),
            confirmed_or_possible_rate=("confirmation_binary_moderate", "mean"),
            mean_support_score=("anomaly_support_score", "mean"),
        )
        .reset_index()
    )
    anomaly_confidence_summary["confirmed_or_possible_percent"] = (
        anomaly_confidence_summary["confirmed_or_possible_rate"] * 100
    )
    outputs["anomaly_confidence_summary"] = anomaly_confidence_summary

    # Direction confidence vs rating-level direction agreement.
    direction_confidence = merged[merged["direction_match_binary"].notna()].copy()
    if not direction_confidence.empty:
        direction_confidence["direction_match_binary"] = direction_confidence[
            "direction_match_binary"
        ].astype(bool)
        outputs["direction_confidence_vs_match"] = summarize_binary_rate(
            direction_confidence,
            "direction_match_binary",
            group_cols=["direction_confidence_normalized"],
            label="direction_agreement_rate",
        )

    return outputs


# =============================================================================
# RTK-GUIDED EXTENDED ANALYSES USING JOINED ANOMALY ATTRIBUTES
# =============================================================================


def summarize_numeric_by_group(
    data: pd.DataFrame,
    group_col: str,
    numeric_cols: list[str],
) -> pd.DataFrame:
    """Summarize numeric anomaly metrics by a categorical group."""
    rows = []
    if data.empty or group_col not in data.columns:
        return pd.DataFrame()
    for group_value, group in data.groupby(group_col, dropna=False):
        for col in numeric_cols:
            if col not in group.columns:
                continue
            numeric = pd.to_numeric(group[col], errors="coerce").dropna()
            if numeric.empty:
                continue
            rows.append(
                {
                    group_col: group_value,
                    "metric": col,
                    "n": int(numeric.size),
                    "mean": float(numeric.mean()),
                    "median": float(numeric.median()),
                    "std": float(numeric.std()) if numeric.size > 1 else np.nan,
                    "min": float(numeric.min()),
                    "max": float(numeric.max()),
                }
            )
    return pd.DataFrame(rows)


def add_quantile_bin(data: pd.DataFrame, source_col: str, output_col: str, q: int = 3) -> pd.DataFrame:
    """Add a quantile-bin column when enough unique numeric values exist."""
    out = data.copy()
    if source_col not in out.columns:
        return out
    numeric = pd.to_numeric(out[source_col], errors="coerce")
    if numeric.notna().sum() < q or numeric.nunique(dropna=True) < q:
        out[output_col] = pd.NA
        return out
    try:
        out[output_col] = pd.qcut(numeric, q=q, duplicates="drop")
    except Exception:
        out[output_col] = pd.NA
    return out


def run_rtk_extended_anomaly_attribute_analysis(
    df_rtk: pd.DataFrame,
    consensus_rtk: pd.DataFrame,
    pairwise_evidence: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Additional RTK-guided analyses enabled by the joined anomaly attributes.

    These analyses are intended as manuscript-supporting diagnostics, not as a
    replacement for the core blind/RTK confirmation-rate and direction-agreement
    summaries. They use z-score for anomaly direction/severity, while retaining
    anomaly area, crop, and issue categories as contextual or extent variables.
    Residuals are intentionally excluded from RTK direction and severity outputs.
    """
    outputs: dict[str, pd.DataFrame] = {}
    if df_rtk.empty or consensus_rtk.empty:
        return outputs

    moderate_col = "moderate_confirmed_3of4_or_75pct_confirmed_possible"
    rtk = consensus_rtk.copy()
    if moderate_col in rtk.columns:
        rtk["rtk_confirmed_moderate"] = rtk[moderate_col].astype("boolean")

    # Basic joined-anomaly and pairwise diagnostics.
    diagnostic_rows = [
        {
            "metric": "rtk_analysis_rows",
            "value": len(df_rtk),
        },
        {
            "metric": "rtk_point_level_evidence_rows",
            "value": len(consensus_rtk),
        },
        {
            "metric": "rtk_pairwise_evidence_rows",
            "value": len(pairwise_evidence) if isinstance(pairwise_evidence, pd.DataFrame) else 0,
        },
        {
            "metric": "rtk_rows_with_joined_anomaly_id",
            "value": int(df_rtk.get("rtk_joined_anomaly_ids", pd.Series(dtype=object)).astype(str).str.len().gt(0).sum()) if "rtk_joined_anomaly_ids" in df_rtk.columns else np.nan,
        },
        {
            "metric": "rtk_rows_with_evaluable_direction",
            "value": int(df_rtk.get("direction_match_binary", pd.Series(dtype=object)).notna().sum()) if "direction_match_binary" in df_rtk.columns else np.nan,
        },
    ]
    outputs["rtk_joined_anomaly_diagnostics"] = pd.DataFrame(diagnostic_rows)

    # Confusion-style table: observed RTK direction vs STADS direction.
    if {"consensus_direction_mode", "stads_direction_normalized"}.issubset(rtk.columns):
        confusion_df = rtk[
            rtk["consensus_direction_mode"].isin(["Positive", "Negative", "Same"])
            & rtk["stads_direction_normalized"].isin(["Positive", "Negative"])
        ].copy()
        if not confusion_df.empty:
            outputs["rtk_direction_confusion"] = pd.crosstab(
                confusion_df["stads_direction_normalized"],
                confusion_df["consensus_direction_mode"],
                margins=True,
            ).reset_index()

    # Confirmation by issue categories noted in the field.
    issue_cols = [col for col in RTK_ISSUE_COLS if col in rtk.columns]
    if issue_cols and moderate_col in rtk.columns:
        issue_long = []
        for col in issue_cols:
            tmp = rtk[[col, moderate_col]].copy()
            tmp = tmp.rename(columns={col: "issue_category"})
            tmp["issue_source_col"] = col
            tmp = tmp.dropna(subset=["issue_category"])
            if not tmp.empty:
                issue_long.append(tmp)
        if issue_long:
            issues = pd.concat(issue_long, ignore_index=True)
            issues["issue_category"] = issues["issue_category"].astype(str).str.strip()
            outputs["rtk_confirmation_by_issue"] = summarize_binary_rate(
                issues,
                moderate_col,
                group_cols=["issue_source_col", "issue_category"],
                label="confirmation_rate",
            )

    # Confirmation and direction agreement by anomaly metric bins.
    metric_bin_specs = [
        ("abs_stads_zscore", "abs_zscore_bin"),
        ("anomaly_area_m2", "anomaly_area_bin"),
    ]
    for metric_col, bin_col in metric_bin_specs:
        if metric_col not in rtk.columns:
            continue
        binned = add_quantile_bin(rtk, metric_col, bin_col, q=3)
        if bin_col not in binned.columns or binned[bin_col].isna().all():
            continue
        binned[bin_col] = binned[bin_col].astype(str)
        if moderate_col in binned.columns:
            outputs[f"rtk_confirmation_by_{bin_col}"] = summarize_binary_rate(
                binned.dropna(subset=[moderate_col]),
                moderate_col,
                group_cols=[bin_col],
                label="confirmation_rate",
            )
        if "direction_consensus_match" in binned.columns:
            dir_eval = binned[binned["direction_consensus_match"].notna()].copy()
            if not dir_eval.empty:
                dir_eval["direction_consensus_match"] = dir_eval["direction_consensus_match"].astype(bool)
                outputs[f"rtk_direction_agreement_by_{bin_col}"] = summarize_binary_rate(
                    dir_eval,
                    "direction_consensus_match",
                    group_cols=[bin_col],
                    label="direction_agreement_rate",
                )

    # Numeric anomaly metric summaries by confirmation and direction-match status.
    numeric_cols = [
        "stads_zscore",
        "abs_stads_zscore",
                "anomaly_area_m2",
        "anomaly_perimeter_m",
        "ndvi",
        "Predicted_NDVI",
        "cum_gdd",
        "cum_precip",
    ]
    numeric_cols = [col for col in numeric_cols if col in rtk.columns]
    if numeric_cols and moderate_col in rtk.columns:
        rtk["rtk_confirmation_status_moderate"] = np.where(
            rtk[moderate_col].astype(bool),
            "Confirmed_or_possible",
            "Not_confirmed",
        )
        outputs["rtk_anomaly_metrics_by_confirmation"] = summarize_numeric_by_group(
            rtk,
            group_col="rtk_confirmation_status_moderate",
            numeric_cols=numeric_cols,
        )
    if numeric_cols and "direction_consensus_match" in rtk.columns:
        dir_metric = rtk[rtk["direction_consensus_match"].notna()].copy()
        if not dir_metric.empty:
            dir_metric["rtk_direction_match_status"] = np.where(
                dir_metric["direction_consensus_match"].astype(bool),
                "Match",
                "Mismatch",
            )
            outputs["rtk_anomaly_metrics_by_direction_match"] = summarize_numeric_by_group(
                dir_metric,
                group_col="rtk_direction_match_status",
                numeric_cols=numeric_cols,
            )

    # Joined anomaly attribute completeness at the cleaned RTK row level.
    completeness_cols = [
        "OBJECTID",
        "point_id",
        "stads_zscore",
            "anomaly_area_m2",
        "CropName",
        "imagerydate",
        "bad_things",
        "bad_thin_1",
    ]
    completeness_cols = [col for col in completeness_cols if col in df_rtk.columns]
    if completeness_cols:
        outputs["rtk_joined_attribute_completeness"] = (
            df_rtk[completeness_cols]
            .notna()
            .sum()
            .reset_index()
            .rename(columns={"index": "column", 0: "non_null_rows"})
        )
        outputs["rtk_joined_attribute_completeness"]["total_rtk_rows"] = len(df_rtk)
        outputs["rtk_joined_attribute_completeness"]["non_null_percent"] = (
            outputs["rtk_joined_attribute_completeness"]["non_null_rows"] / len(df_rtk) * 100
        )

    return outputs


# =============================================================================
# PLOTTING
# =============================================================================


def save_current_plot(filename: str) -> None:
    """Save the active matplotlib figure as PNG and SVG."""
    if SAVE_FIGURES:
        plt.savefig(OUTPUT_DIR / f"{filename}.png", dpi=FIGURE_DPI, bbox_inches="tight")
        plt.savefig(OUTPUT_DIR / f"{filename}.svg", bbox_inches="tight")
    plt.close()


def plot_stacked_by_crop(
    data: pd.DataFrame,
    category_col: str,
    crop_col: str = "Crop",
    normalize: bool = False,
    category_order: Optional[list[str]] = None,
    crop_order: Optional[list[str]] = None,
    title: Optional[str] = None,
    ylabel: Optional[str] = None,
    filename: Optional[str] = None,
) -> pd.DataFrame:
    """Create stacked bar plot of a categorical classification variable by crop."""
    plot_df = data[[crop_col, category_col]].dropna().copy()

    counts = plot_df.groupby([crop_col, category_col]).size().unstack(fill_value=0)

    if crop_order is not None:
        counts = counts.reindex(crop_order)

    if category_order is not None:
        existing_categories = [c for c in category_order if c in counts.columns]
        remaining_categories = [c for c in counts.columns if c not in existing_categories]
        counts = counts[existing_categories + remaining_categories]

    values = counts.div(counts.sum(axis=1), axis=0) * 100 if normalize else counts

    ax = values.plot(
        kind="bar",
        stacked=True,
        figsize=(7.0, 4.5),
        width=0.75,
        edgecolor="black",
        linewidth=0.4,
    )
    ax.set_title(title if title else category_col, fontsize=12)
    ax.set_xlabel("Crop", fontsize=11)
    ax.set_ylabel(
        ylabel if ylabel else ("Percent of observations" if normalize else "Number of observations"),
        fontsize=11,
    )
    ax.tick_params(axis="x", labelrotation=0)
    ax.tick_params(axis="both", labelsize=10)
    if normalize:
        ax.set_ylim(0, 100)
    ax.legend(
        title=category_col,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=9,
        title_fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    if filename:
        save_current_plot(filename)

    return counts


def plot_rate_with_ci(
    summary: pd.DataFrame,
    x_col: str,
    rate_col: str,
    title: str,
    ylabel: str,
    filename: str,
    hue_col: Optional[str] = None,
) -> None:
    """Simple point/bar plot for rates with Wilson confidence intervals."""
    if summary.empty or rate_col not in summary.columns:
        return

    plot_df = summary.copy()
    plot_df = plot_df.dropna(subset=[rate_col])
    if plot_df.empty:
        return

    if hue_col and hue_col in plot_df.columns:
        # Keep this simple: one line/bar group per hue using dodged x positions.
        labels = plot_df[x_col].astype(str).tolist()
        unique_hues = plot_df[hue_col].astype(str).unique().tolist()
        fig, ax = plt.subplots(figsize=(max(7.0, len(labels) * 0.45), 4.5))
        x_base = np.arange(plot_df[x_col].nunique())
        width = 0.8 / max(1, len(unique_hues))

        x_labels = plot_df[x_col].drop_duplicates().astype(str).tolist()
        x_lookup = {label: i for i, label in enumerate(x_labels)}

        for hue_i, hue in enumerate(unique_hues):
            sub = plot_df[plot_df[hue_col].astype(str) == hue].copy()
            x = np.array([x_lookup[str(v)] for v in sub[x_col]]) + (hue_i - (len(unique_hues) - 1) / 2) * width
            y = sub[rate_col].astype(float) * 100
            yerr = np.vstack(
                [
                    y - sub["ci_lower"].astype(float) * 100,
                    sub["ci_upper"].astype(float) * 100 - y,
                ]
            )
            ax.bar(x, y, width=width, edgecolor="black", linewidth=0.4, label=hue)
            ax.errorbar(x, y, yerr=yerr, fmt="none", capsize=3, linewidth=1)

        ax.set_xticks(x_base)
        ax.set_xticklabels(x_labels, rotation=30, ha="right")
        ax.legend(title=hue_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    else:
        fig, ax = plt.subplots(figsize=(max(7.0, len(plot_df) * 0.55), 4.5))
        x = np.arange(len(plot_df))
        y = plot_df[rate_col].astype(float) * 100
        yerr = np.vstack(
            [
                y - plot_df["ci_lower"].astype(float) * 100,
                plot_df["ci_upper"].astype(float) * 100 - y,
            ]
        )
        ax.bar(x, y, edgecolor="black", linewidth=0.4)
        ax.errorbar(x, y, yerr=yerr, fmt="none", capsize=3, linewidth=1)
        ax.set_xticks(x)
        ax.set_xticklabels(plot_df[x_col].astype(str), rotation=30, ha="right")

    ax.set_title(title, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_ylim(0, 100)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save_current_plot(filename)


def run_plots(df: pd.DataFrame, consensus_outputs: dict[str, pd.DataFrame], direction_outputs: dict[str, pd.DataFrame], confidence_outputs: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Generate existing and new publication-oriented plots."""
    plot_tables: dict[str, pd.DataFrame] = {}

    crop_order = df["Crop"].dropna().value_counts().index.tolist() if "Crop" in df.columns else None

    plot_tables["anomaly_by_crop_counts"] = plot_stacked_by_crop(
        df,
        category_col="anomaly_indicated_normalized",
        crop_order=crop_order,
        category_order=ANOMALY_INDICATED_ORDER,
        normalize=False,
        title="Anomaly Classification by Crop",
        ylabel="Number of reviewer classifications",
        filename="anomaly_indicated_by_crop_counts",
    )
    plot_tables["anomaly_by_crop_percent"] = plot_stacked_by_crop(
        df,
        category_col="anomaly_indicated_normalized",
        crop_order=crop_order,
        category_order=ANOMALY_INDICATED_ORDER,
        normalize=True,
        title="Anomaly Classification by Crop",
        ylabel="Percent of reviewer classifications",
        filename="anomaly_indicated_by_crop_percent",
    )
    plot_tables["direction_match_by_crop_counts"] = plot_stacked_by_crop(
        df,
        category_col="direction_match_status",
        crop_order=crop_order,
        category_order=["Match", "Mismatch", "Indeterminate"],
        normalize=False,
        title="Agreement Between Classified Direction and STADS Direction by Crop",
        ylabel="Number of reviewer classifications",
        filename="direction_match_by_crop_counts",
    )
    plot_tables["direction_match_by_crop_percent"] = plot_stacked_by_crop(
        df,
        category_col="direction_match_status",
        crop_order=crop_order,
        category_order=["Match", "Mismatch", "Indeterminate"],
        normalize=True,
        title="Agreement Between Classified Direction and STADS Direction by Crop",
        ylabel="Percent of reviewer classifications",
        filename="direction_match_by_crop_percent",
    )

    if "consensus_rule_sensitivity" in consensus_outputs:
        plot_rate_with_ci(
            consensus_outputs["consensus_rule_sensitivity"],
            x_col="consensus_rule",
            rate_col="confirmation_rate",
            title="Consensus Confirmation Rate by Rule",
            ylabel="Confirmation rate (%)",
            filename="consensus_confirmation_rate_by_rule",
        )

    if "direction_consensus_by_crop" in direction_outputs:
        plot_rate_with_ci(
            direction_outputs["direction_consensus_by_crop"],
            x_col="Crop",
            rate_col="direction_agreement_rate",
            title="Consensus Direction Agreement by Crop",
            ylabel="Direction agreement (%)",
            filename="consensus_direction_agreement_by_crop",
        )

    if "anomaly_confidence_vs_consensus" in confidence_outputs:
        plot_rate_with_ci(
            confidence_outputs["anomaly_confidence_vs_consensus"],
            x_col="anomaly_confidence_normalized",
            rate_col="consensus_agreement_rate",
            title="Reviewer Confidence vs Agreement with Consensus",
            ylabel="Agreement with consensus (%)",
            filename="anomaly_confidence_vs_consensus_agreement",
        )

    return plot_tables


# =============================================================================
# OUTPUT WRITING
# =============================================================================


def safe_sheet_name(name: str, used: set[str]) -> str:
    """Excel sheet names must be <=31 chars and unique."""
    base = name[:31]
    candidate = base
    i = 1
    while candidate in used:
        suffix = f"_{i}"
        candidate = base[: 31 - len(suffix)] + suffix
        i += 1
    used.add(candidate)
    return candidate


def write_tables_to_excel(tables: dict[str, pd.DataFrame], output_path: Path) -> None:
    """Write a dict of DataFrames to an Excel workbook."""
    used_names: set[str] = set()
    with pd.ExcelWriter(output_path) as writer:
        for name, table in tables.items():
            if isinstance(table, pd.DataFrame):
                sheet = safe_sheet_name(name, used_names)
                table.to_excel(writer, sheet_name=sheet, index=False)


def write_model_outputs(model_outputs: dict[str, object], output_dir: Path) -> None:
    """Write model outputs to text/CSV files."""
    summary_text = model_outputs.get("mixed_model_summary_text")
    if summary_text:
        (output_dir / "mixed_effects_or_clustered_model_summary.txt").write_text(
            str(summary_text),
            encoding="utf-8",
        )

    model_metadata = {
        k: v
        for k, v in model_outputs.items()
        if k != "mixed_model_summary_text"
    }
    if model_metadata:
        pd.DataFrame(
            [{k: str(v) for k, v in model_metadata.items()}]
        ).to_csv(output_dir / "model_metadata.csv", index=False)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df_locations, in_anomaly_lookup = build_location_lookup()

    # -------------------------------------------------------------------------
    # Load blind independent-reviewer evidence and optional RTK-guided evidence.
    # -------------------------------------------------------------------------
    df_blind_raw, df_filtered_out = load_reviewer_spreadsheets(in_anomaly_lookup)
    df_rtk_raw, df_rtk_filtered_out, df_rtk_pairwise_evidence = load_rtk_guided_gpkg(
        df_locations,
        in_anomaly_lookup,
    )

    raw_frames = [df_blind_raw]
    if not df_rtk_raw.empty:
        raw_frames.append(df_rtk_raw)
    df_all_raw = pd.concat(raw_frames, ignore_index=True, sort=False)
    df_all = clean_and_derive_rating_columns(df_all_raw)

    df_blind = df_all[df_all["evaluation_method_analysis"] == BLIND_METHOD_LABEL].copy()
    df_rtk = df_all[df_all["evaluation_method_analysis"] == RTK_METHOD_LABEL].copy()

    if len(df_filtered_out) > 0:
        filtered_out_output_path = OUTPUT_DIR / "blind_classification_rows_filtered_out_not_in_anomaly.xlsx"
        df_filtered_out.to_excel(filtered_out_output_path, index=False)
        print(f"Saved blind filtered-out row audit file to: {filtered_out_output_path}")

    if len(df_rtk_filtered_out) > 0:
        rtk_filtered_out_output_path = OUTPUT_DIR / "rtk_rows_filtered_out_not_in_anomaly.xlsx"
        df_rtk_filtered_out.to_excel(rtk_filtered_out_output_path, index=False)
        print(f"Saved RTK filtered-out row audit file to: {rtk_filtered_out_output_path}")

    if len(df_rtk_pairwise_evidence) > 0:
        rtk_pairwise_output_path = OUTPUT_DIR / "rtk_pairwise_compare_to_inference_audit.xlsx"
        df_rtk_pairwise_evidence.to_excel(rtk_pairwise_output_path, index=False)
        print(f"Saved RTK pairwise compare_to inference audit file to: {rtk_pairwise_output_path}")

    # -------------------------------------------------------------------------
    # Core analyses.
    # Reliability, reviewer confidence, and mixed modeling are run on blind
    # independent-reviewer ratings only. RTK-guided evidence is a single-observer
    # stream and is therefore used for method comparison and direction agreement,
    # not inter-rater reliability.
    # -------------------------------------------------------------------------
    reliability_outputs = run_inter_rater_reliability(df_blind)

    consensus_blind = build_point_level_consensus(df_blind)
    consensus_rtk = build_point_level_consensus(df_rtk) if not df_rtk.empty else pd.DataFrame()
    consensus_all = pd.concat(
        [c for c in [consensus_blind, consensus_rtk] if not c.empty],
        ignore_index=True,
        sort=False,
    )

    # -------------------------------------------------------------------------
    # Parallel blind and RTK STADS-support analyses.
    # The underlying functions return generic table names; rename them here so
    # workbook outputs clearly separate the evidence streams.
    # -------------------------------------------------------------------------
    blind_consensus_outputs_raw = run_consensus_support_analysis(consensus_blind)
    blind_consensus_outputs = rename_consensus_outputs(
        blind_consensus_outputs_raw,
        prefix="blind",
    )

    if not consensus_rtk.empty:
        rtk_consensus_outputs_raw = run_consensus_support_analysis(consensus_rtk)
        rtk_consensus_outputs = rename_consensus_outputs(
            rtk_consensus_outputs_raw,
            prefix="rtk",
        )
    else:
        rtk_consensus_outputs_raw = {}
        rtk_consensus_outputs = {}

    blind_direction_outputs_raw = run_direction_agreement_analysis(df_blind, consensus_blind)
    blind_direction_outputs = rename_direction_outputs(
        blind_direction_outputs_raw,
        prefix="blind",
    )

    if not df_rtk.empty and not consensus_rtk.empty:
        rtk_direction_outputs_raw = run_direction_agreement_analysis(df_rtk, consensus_rtk)
        rtk_direction_outputs = rename_direction_outputs(
            rtk_direction_outputs_raw,
            prefix="rtk",
        )
    else:
        rtk_direction_outputs_raw = {}
        rtk_direction_outputs = {}

    # Additional RTK-only analyses enabled by the joined anomaly attributes
    # (z-score severity/direction, anomaly area, crop, and field issue categories).
    rtk_extended_outputs = run_rtk_extended_anomaly_attribute_analysis(
        df_rtk,
        consensus_rtk,
        df_rtk_pairwise_evidence,
    )

    # Combined method-level outputs are retained separately so the workbook has
    # both parallel method-specific tables and an explicit method comparison.
    method_direction_outputs_raw = run_direction_agreement_analysis(df_all, consensus_all)
    method_direction_outputs = rename_direction_outputs(
        method_direction_outputs_raw,
        prefix="method",
    )

    model_outputs = run_mixed_effects_model(df_blind)
    method_outputs = run_method_comparison(consensus_all)
    confidence_outputs = run_confidence_analysis(df_blind, consensus_blind)

    parallel_summary_index = make_parallel_method_summary_index(
        blind_consensus_outputs=blind_consensus_outputs,
        rtk_consensus_outputs=rtk_consensus_outputs,
        blind_direction_outputs=blind_direction_outputs,
        rtk_direction_outputs=rtk_direction_outputs,
        method_outputs=method_outputs,
    )

    # Plots use blind reviewer ratings for the reviewer-classification figures,
    # preserving comparability with the original script. The unprefixed raw
    # outputs are used internally because run_plots expects generic names.
    plot_tables = run_plots(
        df_blind,
        blind_consensus_outputs_raw,
        blind_direction_outputs_raw,
        confidence_outputs,
    )

    # Existing summary statistics from the original script, updated to normalized labels.
    classification_summary = (
        df_blind.groupby(
            by=[
                "anomaly_indicated_normalized",
                "anomaly_confidence_normalized",
                "Crop",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="n")
    )
    total_n = int(classification_summary["n"].sum())
    print(f"Total blind reviewer classifications in classification_summary: {total_n}")
    print(f"RTK-guided in-anomaly field assessments: {len(df_rtk)}")

    # Assemble output workbook. Put the high-level parallel tables first so the
    # workbook is easier to interpret before the detailed audit/model sheets.
    all_tables: dict[str, pd.DataFrame] = {
        "parallel_output_index": parallel_summary_index,
        "all_cleaned_ratings": df_all,
        "blind_cleaned_ratings": df_blind,
        "blind_point_level_consensus": consensus_blind,
        "all_method_point_evidence": consensus_all,
        "classification_summary": classification_summary,
    }
    if not df_rtk.empty:
        all_tables["rtk_cleaned_ratings"] = df_rtk
    if not consensus_rtk.empty:
        all_tables["rtk_point_level_evidence"] = consensus_rtk
    if not df_rtk_pairwise_evidence.empty:
        all_tables["rtk_pairwise_compare_to_audit"] = df_rtk_pairwise_evidence

    # Parallel STADS-support outputs.
    all_tables.update(blind_consensus_outputs)
    all_tables.update(rtk_consensus_outputs)
    all_tables.update(method_outputs)
    all_tables.update(blind_direction_outputs)
    all_tables.update(rtk_direction_outputs)
    all_tables.update(method_direction_outputs)

    # Reliability, confidence, and figure-source tables.
    all_tables.update(reliability_outputs)
    all_tables.update(confidence_outputs)
    all_tables.update(plot_tables)

    summary_output_path = OUTPUT_DIR / "stads_reviewer_and_rtk_performance_analysis_tables.xlsx"
    write_tables_to_excel(all_tables, summary_output_path)
    write_model_outputs(model_outputs, OUTPUT_DIR)

    # Write a short readme for missing inputs / interpretation.
    readme_lines = [
        "STADS reviewer + RTK-guided performance analysis outputs",
        "======================================================",
        "",
        f"Blind in-anomaly reviewer ratings: {len(df_blind)}",
        f"Blind point-level consensus rows: {len(consensus_blind)}",
        f"RTK-guided in-anomaly field-assessment rows: {len(df_rtk)}",
        f"Combined method-level point/evidence rows: {len(consensus_all)}",
        f"Output workbook: {summary_output_path.name}",
        "",
        "Main parallel STADS-support tables:",
        "- blind_confirmation_rates and rtk_confirmation_rates summarize whether STADS-detected anomalous areas were supported by observations under the same confirmation rules.",
        "- method_confirmation_rates directly compares Blind paired vs RTK guided confirmation rates using the moderate rule.",
        "- blind_direction_agreement and rtk_direction_agreement summarize whether observed direction matched STADS direction.",
        "- *_by_crop and *_by_stads_direction tables provide the corresponding stratified summaries.",
        "- parallel_output_index lists the main parallel output tables and whether each was created in this run.",
        "",
        "Important interpretation notes:",
        "- Inter-rater reliability is calculated only for the blind independent-reviewer spreadsheets.",
        "- RTK-guided scouting is treated as a separate single-observer evidence stream, not as another reviewer.",
        "- The RTK GeoPackage path is configured to use STADS_RTK_Anomalies_All_Merged_epsg4326.gpkg and the 'Anomalies With Observations' layer.",
        "- For RTK rows, in-anomaly status is derived from the joined anomaly attributes; STADS direction is derived only from zscore/stads_zscore.",
        "- For RTK zscore/stads_zscore, negative values indicate positive deviance and positive values indicate negative deviance.",
        "- Residuals are intentionally not used for RTK anomaly direction or severity analyses.",
        "- RTK sample_pai is treated as the point-pair identifier; a_or_b identifies A/B within each pair.",
        "- RTK compare_to values are interpreted within each sample_pai pair: Better=Positive, Worse=Negative, Same=Same for the current A/B point relative to the paired A/B point.",
        "- If compare_to is recorded on the out-of-anomaly paired point, the direction is inverted before assigning evidence to the in-anomaly point.",
        "- If compare_to is recorded on both A and B, both entries are used as pairwise evidence; conflicting inferred directions become Indeterminate/Low confidence after collapse.",
        "- For RTK anomaly presence, inferred Positive/Negative are treated as Confirmed and Same as Not Confirmed.",
        "- Directional correctness is evaluated separately by comparing the mapped RTK direction to STADS direction.",
        "- Krippendorff's alpha / kappa / AC1 assess reviewer consistency, not STADS performance by themselves.",
        "- Consensus confirmation rates and direction agreement are the main STADS-support summaries.",
        "- Additional RTK tables summarize confirmation and direction agreement by abs(zscore), anomaly area, and field issue categories; residual-based severity tables are not generated.",
        "- The mixed model is exploratory unless the candidate covariates and grouping variables are finalized.",
        "",
    ]

    if "method_comparison_note" in method_outputs:
        readme_lines.append(method_outputs["method_comparison_note"].iloc[0]["note"])
        readme_lines.append("")

    if model_outputs.get("model_note"):
        readme_lines.append(str(model_outputs["model_note"]))
        readme_lines.append("")

    if model_outputs.get("mixed_model_error"):
        readme_lines.append("Mixed-model attempt error:")
        readme_lines.append(str(model_outputs["mixed_model_error"]))
        readme_lines.append("")

    (OUTPUT_DIR / "README_analysis_notes.txt").write_text("\n".join(readme_lines), encoding="utf-8")

    print(f"Saved tables, figures, and notes to: {OUTPUT_DIR}")

#%%
if __name__ == "__main__":
    main()

# %%
