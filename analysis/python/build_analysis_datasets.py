from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "output" / "analysis"

GROUP_FILE = ROOT / "分组情况9.20.xlsx"
RCT_FILE = ROOT / "Counsellor RCT.xlsx"
BEHAVIOR_FILE = ROOT / "behavior.xlsx"
PRACTICE_FILE = ROOT / "practice_data_summary.xlsx"


def extract_token(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    text = re.sub(r"\.0$", "", text)
    match = re.search(r"[A-Za-z0-9]+", text)
    return match.group(0) if match else text


def extract_digits(value: object) -> str | None:
    token = extract_token(value)
    if token is None:
        return None
    match = re.search(r"\d+", token)
    return match.group(0) if match else None


def trailing_three(value: object) -> str | None:
    token = extract_token(value)
    if token is None:
        return None
    digits = re.fullmatch(r"\d+", token)
    if digits:
        return digits.group(0)[-3:].zfill(3)
    return token


def score_block(df: pd.DataFrame, columns: Iterable[str], prefix: str) -> pd.DataFrame:
    columns = list(columns)
    if not columns:
        return pd.DataFrame(index=df.index)
    numeric = df[columns].apply(pd.to_numeric, errors="coerce")
    answered = numeric.notna().sum(axis=1)
    totals = numeric.sum(axis=1, min_count=1)
    means = numeric.mean(axis=1)
    return pd.DataFrame(
        {
            f"{prefix}_items_answered": answered,
            f"{prefix}_total": totals,
            f"{prefix}_mean": means,
        },
        index=df.index,
    )


def load_roster() -> tuple[pd.DataFrame, dict[int, str]]:
    roster = pd.read_excel(
        GROUP_FILE,
        usecols=["分组", "我的学号是______（仅用于问卷编号）。", "前测完成情况"],
    ).rename(
        columns={
            "分组": "group_label",
            "我的学号是______（仅用于问卷编号）。": "full_id_raw",
            "前测完成情况": "pretest_complete_status",
        }
    )
    roster["full_id"] = roster["full_id_raw"].map(extract_digits)
    roster["short_id"] = roster["full_id"].str[-3:]

    t1 = pd.read_excel(RCT_FILE, sheet_name="T1", usecols=["ID", "Group"]).rename(
        columns={"Group": "group_code"}
    )
    t1["full_id"] = t1["ID"].map(extract_digits)
    code_map = (
        t1.merge(roster[["full_id", "group_label"]], on="full_id", how="left")
        .dropna(subset=["group_label"])
        .drop_duplicates(subset=["group_code", "group_label"])
    )
    if code_map["group_code"].nunique() != code_map["group_label"].nunique():
        raise ValueError("Group code mapping is not one-to-one.")
    group_code_to_label = dict(
        sorted(code_map[["group_code", "group_label"]].itertuples(index=False, name=None))
    )
    group_label_to_code = {label: code for code, label in group_code_to_label.items()}
    roster["group_code"] = roster["group_label"].map(group_label_to_code)
    return roster, group_code_to_label


def build_rct_outputs(roster: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    scale_prefixes = ["HS", "SM", "CC", "Eff", "Use"]
    sheet_order = ["T1", "T2", "T3", "T4", "T5", "T6"]

    for idx, sheet in enumerate(sheet_order, start=1):
        df = pd.read_excel(RCT_FILE, sheet_name=sheet)
        df["full_id"] = df["ID"].map(extract_digits)
        df["_source_order"] = range(len(df))
        df["submission_time"] = pd.to_datetime(df.get("Time"), errors="coerce")
        df = df.sort_values(["full_id", "submission_time", "_source_order"]).copy()
        if sheet == "T6":
            sav_override_ids = {"<STUDENT_ID_047>"}
            special = df[df["full_id"].isin(sav_override_ids)].drop_duplicates(
                subset=["full_id"], keep="first"
            )
            regular = df[~df["full_id"].isin(sav_override_ids)].drop_duplicates(
                subset=["full_id"], keep="last"
            )
            df = pd.concat([regular, special], ignore_index=True)
            df = df.sort_values(["full_id", "submission_time", "_source_order"]).reset_index(
                drop=True
            )
        else:
            df = df.drop_duplicates(subset=["full_id"], keep="last").reset_index(drop=True)

        out = pd.DataFrame(index=df.index)
        out["timepoint"] = sheet
        out["timepoint_index"] = idx
        out["tp_raw"] = df.get("TP")
        out["time_raw"] = df.get("Time")
        out["submission_time"] = df.get("submission_time")
        out["full_id"] = df["full_id"]
        out = out.merge(
            roster[["full_id", "short_id", "group_code", "group_label"]],
            on="full_id",
            how="left",
        )

        for prefix in scale_prefixes:
            columns = [col for col in df.columns if str(col).startswith(f"{prefix}{sheet}_")]
            out = pd.concat([out, score_block(df, columns, prefix.lower())], axis=1)

        if sheet == "T6":
            pa_columns = [col for col in df.columns if re.fullmatch(r"PA_\d+", str(col))]
            na_columns = [col for col in df.columns if re.fullmatch(r"NA_\d+", str(col))]
            out = pd.concat([out, score_block(df, pa_columns, "pa")], axis=1)
            out = pd.concat([out, score_block(df, na_columns, "na")], axis=1)
            if "Sex" in df.columns:
                out["sex_t6"] = df["Sex"]
            if "Age" in df.columns:
                out["age_t6"] = pd.to_numeric(df["Age"], errors="coerce")

        out["self_efficacy_total"] = out[["hs_total", "sm_total", "cc_total"]].sum(
            axis=1, min_count=1
        )
        out["self_efficacy_items_answered"] = out[
            ["hs_items_answered", "sm_items_answered", "cc_items_answered"]
        ].sum(axis=1, min_count=1)
        out["self_efficacy_mean"] = (
            out["self_efficacy_total"] / out["self_efficacy_items_answered"]
        )
        frames.append(out)

    rct_long = pd.concat(frames, ignore_index=True).sort_values(["timepoint_index", "full_id"])
    wide = rct_long.pivot_table(
        index=["full_id", "short_id", "group_code", "group_label"],
        columns="timepoint",
        values=[
            "hs_total",
            "hs_mean",
            "sm_total",
            "sm_mean",
            "cc_total",
            "cc_mean",
            "self_efficacy_total",
            "self_efficacy_mean",
            "eff_total",
            "eff_mean",
            "use_total",
            "use_mean",
            "pa_total",
            "pa_mean",
            "na_total",
            "na_mean",
        ],
        aggfunc="first",
    )
    wide.columns = [f"{value}_{timepoint.lower()}" for value, timepoint in wide.columns]
    rct_wide = wide.reset_index().sort_values(["group_code", "full_id"])
    return rct_long, rct_wide


def build_behavior_output(roster: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_excel(BEHAVIOR_FILE)
    df["raw_user_id"] = df["user_id"].map(extract_token)
    df["full_id"] = df["user_id"].map(extract_digits)
    df["short_id_from_behavior"] = df["user_id"].map(trailing_three)

    exact = roster[["full_id", "short_id", "group_code", "group_label"]].drop_duplicates()
    out = df.merge(exact, on="full_id", how="left")
    out["roster_match_status"] = out["group_label"].map(
        lambda value: "exact_full_id_match" if pd.notna(value) else "no_roster_match"
    )

    composites = {
        "emotional_support": [
            "Affirmation and Reassurance",
            "Reflection of Feelings",
        ],
        "insight": [
            "Interpretation",
            "Confrontation",
            "Invite to Take New Perspectives",
        ],
        "content_exploration": [
            "Inquiring Objective Information",
            "Inquiring Subjective Information",
            "Restatement",
            "Minimal Encouragement",
        ],
        "action": [
            "Answer",
            "Invite to Explore New Actions",
        ],
    }

    for metric_type in ["raw", "norm"]:
        for phase in ["pre", "post"]:
            for name, parts in composites.items():
                source_cols = [f"{part}_{metric_type}_{phase}" for part in parts]
                out[f"{name}_{metric_type}_{phase}"] = out[source_cols].apply(
                    pd.to_numeric, errors="coerce"
                ).sum(axis=1, min_count=1)
        for name in composites:
            out[f"{name}_{metric_type}_change"] = (
                out[f"{name}_{metric_type}_post"] - out[f"{name}_{metric_type}_pre"]
            )

    return out


def build_practice_lookup(roster: pd.DataFrame) -> pd.DataFrame:
    unique_ids = sorted(
        set(
            pd.read_excel(PRACTICE_FILE, sheet_name="student_summary", usecols=["student_id"])[
                "student_id"
            ]
            .map(extract_token)
            .dropna()
            .tolist()
        )
    )
    short_counts = roster["short_id"].value_counts(dropna=False)
    rows = []
    for raw_id in unique_ids:
        derived_short_id = trailing_three(raw_id)
        candidates = roster[roster["short_id"] == derived_short_id][
            ["full_id", "short_id", "group_code", "group_label"]
        ].drop_duplicates()

        row = {
            "practice_student_id": raw_id,
            "derived_short_id": derived_short_id,
            "candidate_count": len(candidates),
            "matched_full_id": None,
            "matched_group_code": None,
            "matched_group_label": None,
            "candidate_full_ids": None,
            "candidate_group_labels": None,
        }

        if len(candidates) == 1:
            candidate = candidates.iloc[0]
            row["matched_full_id"] = candidate["full_id"]
            row["matched_group_code"] = candidate["group_code"]
            row["matched_group_label"] = candidate["group_label"]
            row["roster_match_status"] = "unique_short_id_match"
        elif len(candidates) > 1:
            row["candidate_full_ids"] = ";".join(candidates["full_id"].astype(str))
            row["candidate_group_labels"] = ";".join(candidates["group_label"].astype(str))
            row["roster_match_status"] = "ambiguous_short_id_match"
        else:
            row["roster_match_status"] = "no_roster_match"

        rows.append(row)

    return pd.DataFrame(rows).sort_values("practice_student_id")


def build_practice_outputs(
    practice_lookup: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    session_level = pd.read_excel(PRACTICE_FILE, sheet_name="session_level")
    student_summary = pd.read_excel(PRACTICE_FILE, sheet_name="student_summary")

    for frame in [session_level, student_summary]:
        frame["practice_student_id"] = frame["student_id"].map(extract_token)

    session_out = session_level.merge(practice_lookup, on="practice_student_id", how="left")
    student_out = student_summary.merge(practice_lookup, on="practice_student_id", how="left")

    return session_out, student_out


def build_participant_master(
    roster: pd.DataFrame,
    rct_long: pd.DataFrame,
    behavior: pd.DataFrame,
    practice_lookup: pd.DataFrame,
    practice_student: pd.DataFrame,
) -> pd.DataFrame:
    master = roster[
        ["full_id", "short_id", "group_code", "group_label", "pretest_complete_status"]
    ].drop_duplicates()

    rct_presence = (
        rct_long.assign(present=1)
        .pivot_table(
            index="full_id",
            columns="timepoint",
            values="present",
            aggfunc="max",
            fill_value=0,
        )
        .rename(columns=lambda value: f"has_{str(value).lower()}")
        .reset_index()
    )

    rct_selected = rct_long[
        [
            "full_id",
            "timepoint",
            "hs_total",
            "sm_total",
            "cc_total",
            "self_efficacy_total",
            "eff_total",
            "use_total",
            "pa_total",
            "na_total",
            "sex_t6",
            "age_t6",
        ]
    ].copy()
    rct_selected["timepoint"] = rct_selected["timepoint"].str.lower()
    rct_wide = rct_selected.pivot_table(
        index="full_id",
        columns="timepoint",
        values=[
            "hs_total",
            "sm_total",
            "cc_total",
            "self_efficacy_total",
            "eff_total",
            "use_total",
            "pa_total",
            "na_total",
            "sex_t6",
            "age_t6",
        ],
        aggfunc="first",
    )
    rct_wide.columns = [f"{value}_{timepoint}" for value, timepoint in rct_wide.columns]
    rct_wide = rct_wide.reset_index()

    behavior_selected = behavior[
        [
            "full_id",
            "roster_match_status",
            "emotional_support_raw_pre",
            "emotional_support_raw_post",
            "emotional_support_raw_change",
            "insight_raw_pre",
            "insight_raw_post",
            "insight_raw_change",
            "content_exploration_raw_pre",
            "content_exploration_raw_post",
            "content_exploration_raw_change",
            "action_raw_pre",
            "action_raw_post",
            "action_raw_change",
            "emotional_support_norm_pre",
            "emotional_support_norm_post",
            "emotional_support_norm_change",
            "insight_norm_pre",
            "insight_norm_post",
            "insight_norm_change",
            "content_exploration_norm_pre",
            "content_exploration_norm_post",
            "content_exploration_norm_change",
            "action_norm_pre",
            "action_norm_post",
            "action_norm_change",
        ]
    ].rename(columns={"roster_match_status": "behavior_match_status"})

    practice_unique = practice_student[
        practice_student["roster_match_status"] == "unique_short_id_match"
    ].copy()
    practice_unique = practice_unique[
        [
            "matched_full_id",
            "session_count",
            "feedback_count",
            "total_score",
            "valid_scores_n",
            "valid_scores_sum",
            "overall_mean_score_weighted",
            "counselor_count",
        ]
    ].rename(
        columns={
            "matched_full_id": "full_id",
            "session_count": "practice_session_count",
            "feedback_count": "practice_feedback_count",
            "total_score": "practice_total_score",
            "valid_scores_n": "practice_valid_scores_n",
            "valid_scores_sum": "practice_valid_scores_sum",
            "overall_mean_score_weighted": "practice_overall_mean_score_weighted",
            "counselor_count": "practice_counselor_count",
        }
    )

    ambiguous_short_ids = set(
        practice_lookup.loc[
            practice_lookup["roster_match_status"] == "ambiguous_short_id_match",
            "derived_short_id",
        ]
    )
    unique_short_ids = set(
        practice_lookup.loc[
            practice_lookup["roster_match_status"] == "unique_short_id_match",
            "derived_short_id",
        ]
    )

    master = master.merge(rct_presence, on="full_id", how="left")
    master = master.merge(rct_wide, on="full_id", how="left")
    master = master.merge(behavior_selected, on="full_id", how="left")
    master = master.merge(practice_unique, on="full_id", how="left")

    has_cols = [col for col in master.columns if col.startswith("has_t")]
    if has_cols:
        master[has_cols] = master[has_cols].fillna(0).astype(int)

    master["behavior_present"] = master["behavior_match_status"].eq("exact_full_id_match")
    master["practice_unique_match_present"] = master["short_id"].isin(unique_short_ids)
    master["practice_ambiguous_short_id_present"] = master["short_id"].isin(ambiguous_short_ids)

    return master.sort_values(["group_code", "full_id"])


def write_summary(
    roster: pd.DataFrame,
    group_code_to_label: dict[int, str],
    rct_long: pd.DataFrame,
    behavior: pd.DataFrame,
    practice_lookup: pd.DataFrame,
    practice_student: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append("# Analysis Dataset Notes")
    lines.append("")
    lines.append("## Group Mapping")
    lines.append("")
    for code, label in group_code_to_label.items():
        lines.append(f"- Group {code} = {label}")
    lines.append("")
    lines.append("## Roster")
    lines.append("")
    lines.append(f"- Official roster participants: {roster['full_id'].nunique()}")
    duplicate_short = roster[roster["short_id"].duplicated(keep=False)][
        ["full_id", "group_label", "short_id"]
    ].drop_duplicates()
    if duplicate_short.empty:
        lines.append("- Duplicate short IDs in roster: none")
    else:
        dup_text = ", ".join(
            f"{row.short_id}: {row.full_id} ({row.group_label})"
            for row in duplicate_short.itertuples(index=False)
        )
        lines.append(f"- Duplicate short IDs in roster: {dup_text}")
    lines.append("")
    lines.append("## RCT Coverage")
    lines.append("")
    for timepoint in ["T1", "T2", "T3", "T4", "T5", "T6"]:
        n_rows = int((rct_long["timepoint"] == timepoint).sum())
        n_ids = int(rct_long.loc[rct_long["timepoint"] == timepoint, "full_id"].nunique())
        lines.append(f"- {timepoint}: {n_rows} rows, {n_ids} unique IDs")
    raw_duplicates: list[str] = []
    for timepoint in ["T1", "T2", "T3", "T4", "T5", "T6"]:
        raw_df = pd.read_excel(RCT_FILE, sheet_name=timepoint, usecols=["ID"])
        raw_ids = raw_df["ID"].map(extract_digits)
        dup_ids = sorted(raw_ids[raw_ids.duplicated(keep=False)].dropna().unique().tolist())
        if dup_ids:
            raw_duplicates.append(f"{timepoint}: {', '.join(dup_ids)}")
    if raw_duplicates:
        lines.append(
            "- Source duplicate submissions were deduplicated by keeping the latest timestamp, except "
            "T6 ID <STUDENT_ID_047> which keeps the earlier SAV-matched submission: "
            + "; ".join(raw_duplicates)
        )
    lines.append("")
    lines.append("## Behavior Coverage")
    lines.append("")
    lines.append(
        f"- Exact roster-matched behavior rows: "
        f"{int((behavior['roster_match_status'] == 'exact_full_id_match').sum())}"
    )
    unmatched_behavior = sorted(
        behavior.loc[behavior["roster_match_status"] == "no_roster_match", "raw_user_id"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    lines.append(
        "- Behavior IDs with no roster match: "
        + (", ".join(unmatched_behavior) if unmatched_behavior else "none")
    )
    lines.append("")
    lines.append("## Practice Coverage")
    lines.append("")
    match_counts = (
        practice_lookup["roster_match_status"].value_counts(dropna=False).sort_index().to_dict()
    )
    for status, count in match_counts.items():
        lines.append(f"- {status}: {count}")
    unmatched_practice = sorted(
        practice_lookup.loc[
            practice_lookup["roster_match_status"] == "no_roster_match",
            "practice_student_id",
        ]
        .dropna()
        .astype(str)
        .tolist()
    )
    if unmatched_practice:
        lines.append("- Practice IDs with no roster match: " + ", ".join(unmatched_practice))
    ambiguous_practice = practice_lookup[
        practice_lookup["roster_match_status"] == "ambiguous_short_id_match"
    ][["practice_student_id", "candidate_full_ids", "candidate_group_labels"]]
    if not ambiguous_practice.empty:
        for row in ambiguous_practice.itertuples(index=False):
            lines.append(
                f"- Ambiguous practice ID {row.practice_student_id}: "
                f"{row.candidate_full_ids} / {row.candidate_group_labels}"
            )
    unique_practice = practice_student[
        practice_student["roster_match_status"] == "unique_short_id_match"
    ]
    if not unique_practice.empty:
        group_counts = (
            unique_practice["matched_group_label"].value_counts(dropna=False).sort_index().to_dict()
        )
        group_bits = ", ".join(f"{group}: {count}" for group, count in group_counts.items())
        lines.append(f"- Uniquely matched practice students by group: {group_bits}")
    lines.append("")
    lines.append("## Assumptions")
    lines.append("")
    lines.append("- No reverse scoring was applied.")
    lines.append("- HS, SM, CC, Eff, Use, PA, and NA totals are simple item sums.")
    lines.append(
        "- Behavior composites use the user-provided mapping; content exploration includes both objective and subjective information questions."
    )
    lines.append(
        "- Practice files were left conservative: only unique short-ID matches receive participant-level practice metrics."
    )
    (OUTPUT_DIR / "analysis_notes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    roster, group_code_to_label = load_roster()
    rct_long, rct_wide = build_rct_outputs(roster)
    behavior = build_behavior_output(roster)
    practice_lookup = build_practice_lookup(roster)
    practice_session, practice_student = build_practice_outputs(practice_lookup)
    participant_master = build_participant_master(
        roster=roster,
        rct_long=rct_long,
        behavior=behavior,
        practice_lookup=practice_lookup,
        practice_student=practice_student,
    )

    roster.sort_values(["group_code", "full_id"]).to_csv(
        OUTPUT_DIR / "participant_roster.csv", index=False
    )
    rct_long.to_csv(OUTPUT_DIR / "rct_scale_totals_long.csv", index=False)
    rct_wide.to_csv(OUTPUT_DIR / "rct_scale_totals_wide.csv", index=False)
    behavior.sort_values(["group_code", "full_id", "raw_user_id"]).to_csv(
        OUTPUT_DIR / "behavior_clean.csv", index=False
    )
    practice_lookup.to_csv(OUTPUT_DIR / "practice_id_lookup.csv", index=False)
    practice_session.to_csv(OUTPUT_DIR / "practice_session_level_clean.csv", index=False)
    practice_student.to_csv(OUTPUT_DIR / "practice_student_summary_clean.csv", index=False)
    participant_master.to_csv(OUTPUT_DIR / "participant_master.csv", index=False)

    write_summary(
        roster=roster,
        group_code_to_label=group_code_to_label,
        rct_long=rct_long,
        behavior=behavior,
        practice_lookup=practice_lookup,
        practice_student=practice_student,
    )


if __name__ == "__main__":
    main()
