"""Pure functions extracted from the supplied v5 analysis; no raw-data paths."""
from __future__ import annotations
import pandas as pd
import numpy as np
import itertools
from typing import Optional, Iterable

def mode_or_tie(values: Iterable) -> object:
    """Return the unique mode or 'Tie' if multiple values share the top count."""
    cleaned = pd.Series(list(values)).dropna()
    if cleaned.empty:
        return pd.NA
    counts = cleaned.value_counts()
    top = counts[counts == counts.max()].index.tolist()
    return top[0] if len(top) == 1 else "Tie"

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

